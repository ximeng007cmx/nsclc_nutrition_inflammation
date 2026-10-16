# -*- coding: utf-8 -*-
"""make_article_tables.py — 把 results/tables/*.csv 合成文章级表格文档（Markdown + HTML）。

输出：
  results/tables/PathA_文章级表格.md
  results/tables/PathA_文章级表格.html
"""
import csv
import html
import os

TDIR = "D:/WorkBuddy_Y02/out_ali_directions/pathA/results/tables"

TABLES = [
    ("Table1_baseline.csv", "Table 1. Baseline characteristics of the study cohorts",
     "TCGA-LUAD (training cohort) and TCGA-LUSC (cross-histology validation cohort). "
     "Clinical data from TCGA-CDR (Liu J et al., Cell 2018; PMID 29625055). "
     "Signature scores are ssGSEA values of the three nutrition–inflammation axes."),
    ("Table2_axis_immune_correlation.csv",
     "Table 2. Correlation between nutrition–inflammation axes and immune features",
     "Spearman correlations across 1,011 primary tumour samples. "
     "Immune features are MSigDB pathway-level activity proxies "
     "(GO Biological Process and Hallmark; see Table S4b for set sizes and checksums), "
     "NOT deconvolution estimates. "
     "'No overlapping genes' indicates the axis and the immune feature share no genes, "
     "so the correlation is not definitional. FDR = Benjamini–Hochberg."),
    ("Table3_axis_checkpoint_correlation.csv",
     "Table 3. Correlation between nutrition–inflammation axes and immune checkpoint / effector genes",
     "Spearman correlations with single-gene expression (log2(FPKM+1)). "
     "Checkpoint panel: PDCD1, CD274, PDCD1LG2, CTLA4, LAG3, HAVCR2, TIGIT, IDO1, SIGLEC15. "
     "Effector panel: IFNG, GZMA, GZMB, PRF1, CXCL9, CXCL10, CXCL11, HLA-DRA, STAT1, CCR5."),
    ("Table4_internal_validation.csv",
     "Table 4. Internal validation of the prognostic score (negative result)",
     "LUAD split 7:3 at RANDOM (sample.int(); NOT stratified by event status — "
     "observed event rates 38.3% training vs 32.0% held-out, against 36.4% expected "
     "under stratification); feature selection (univariate Cox P<0.10) "
     "and standardisation performed within the training set only. "
     "The held-out C-index is at chance level, indicating the three a-priori axes carry "
     "no prognostic discrimination in LUAD; the previously reported C-index of 0.541 "
     "was an over-optimistic value from fitting on the full dataset."),
    ("Table5_immune_diff_by_nutrition_axis.csv",
     "Table 5. Immune features by nutrition-axis group (independent of definitional overlap)",
     "Samples dichotomised at the median of the nutrition_liver axis, which shares "
     "NO genes with any immune feature. Higher nutrition-axis expression is associated "
     "with higher activity of myeloid/monocyte/macrophage and inflammatory pathways. "
     "Wilcoxon rank-sum test; FDR = Benjamini–Hochberg."),
    ("TableS1_endpoint_sensitivity.csv",
     "Table S1. Sensitivity analysis across endpoints (OS / DSS / PFI)",
     "Same pipeline re-run with each endpoint. LUAD is essentially flat across endpoints "
     "(C-index 0.535–0.542); LUSC shows a significant association for DSS. "
     "Endpoints therefore do not rescue the LUAD training signal."),
    ("TableS2_geneset_sensitivity.csv",
     "Table S2. Sensitivity analysis across feature schemes and training directions",
     "S0 = three a-priori axes (16 genes each); S1 = six MSigDB pathways (1,325 genes); "
     "S2 = 48 prior genes entered directly into LASSO-Cox. "
     "Internal C-index rises with feature number while external C-index falls — "
     "a direct demonstration of overfitting. Cross-histology validation ceiling ≈ 0.55."),
    ("TableS3_immune_diff_by_composite_score.csv",
     "Table S3. Immune features by composite risk score (definitional consistency check)",
     "CAUTION: the composite score contains the lymphocyte axis with a negative coefficient, "
     "so lower immune activity in the high-score group is partly definitional. "
     "This table is provided for completeness only; Table 5 is the non-circular analysis."),
    ("TableS4a_axis_genes.csv",
     "Table S4a. Gene composition of the three a-priori nutrition–inflammation axes",
     "16 genes per axis, entered as a-priori definition (not data-driven selection). "
     "All 48 symbols were detected in the expression matrix. "
     "'detected_in_expression_matrix = N' would indicate a symbol dropped before scoring; "
     "there are none in this dataset."),
    ("TableS4b_immune_features.csv",
     "Table S4b. MSigDB immune features used as pathway-level activity proxies",
     "13 features (10 GO:BP + 3 Hallmark). 'set_size_in_gmt' is the raw MSigDB entry size; "
     "'set_size_in_matrix' is the size after intersecting with the expression matrix. "
     "Overlap columns give the number of genes shared with each axis: the nutrition_liver "
     "axis shares 0 genes with all 13 features, which is what makes its correlations "
     "non-definitional. Source release is MSigDB v2024.1.Hs, established by byte-level "
     "matching against the public MSigDB release archive: md5 of h.all.symbols.gmt "
     "(d3606c325b9862bf3ea7991c2c2bf7d7) and of c5.go.bp.symbols.gmt "
     "(fd535c44825f9be3ea7c5d31e3d45f1d) each match v2024.1.Hs and no other release "
     "between v7.5.1 and v2026.1.Hs."),
    ("TableS4c_checkpoint_effector.csv",
     "Table S4c. Immune checkpoint and effector transcript panel",
     "9 immune checkpoint genes and 10 effector genes, all detected in the expression matrix. "
     "Correlations for these transcripts are reported in Table 3 and Figure 2b."),
    ("TableS4d_msigdb6_pathways.csv",
     "Table S4d. MSigDB pathways used in the S1 feature scheme",
     "6 real MSigDB pathways (source release v2024.1.Hs) spanning 14–732 genes after "
     "intersecting with the expression "
     "matrix (1,325 unique gene-pathway memberships in total). This heterogeneity in set size "
     "is why pathway scores were standardised within the training set before LASSO-Cox."),
]


def read_csv(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.reader(f))


def md_table(rows):
    if not rows:
        return "_(empty)_\n"
    w = len(rows[0])
    out = ["| " + " | ".join(rows[0]) + " |",
           "|" + "|".join(["---"] * w) + "|"]
    for r in rows[1:]:
        r = (r + [""] * w)[:w]
        out.append("| " + " | ".join(x.replace("|", "\\|") for x in r) + " |")
    return "\n".join(out) + "\n"


def html_table(rows):
    if not rows:
        return "<p><em>(empty)</em></p>\n"
    o = ['<table>']
    o.append("<thead><tr>" + "".join("<th>%s</th>" % html.escape(c) for c in rows[0]) + "</tr></thead>")
    o.append("<tbody>")
    for r in rows[1:]:
        o.append("<tr>" + "".join("<td>%s</td>" % html.escape(c) for c in r) + "</tr>")
    o.append("</tbody></table>")
    return "\n".join(o) + "\n"


md = ["# Path A — 文章级表格（营养-炎症分子特征 × 非小细胞肺癌免疫微环境）", "",
      "> 生成日期：2026-10-16  ",
      "> 数据来源：UCSC Xena Toil（TCGA RSEM FPKM，经 log2(FPKM+1) 转换）+ TCGA-CDR 临床表 + MSigDB v2024.1.Hs  ",
      "> 分析样本：TCGA-LUAD 500 例（182 事件）、TCGA-LUSC 491 例（211 事件），均为原发肿瘤  ",
      "> 表述约束：本方案为 ALI 三组分的**分子代理**，不得表述为「等价于 ALI」；",
      "> 免疫特征为**通路级活性代理**，不得表述为 CIBERSORT/ESTIMATE 去卷积结果。", ""]
ht = ["<!DOCTYPE html><html lang='zh'><head><meta charset='utf-8'>",
      "<title>Path A 文章级表格</title><style>",
      "body{font-family:'Segoe UI',system-ui,-apple-system,'Microsoft YaHei',sans-serif;",
      "max-width:1180px;margin:36px auto;padding:0 22px;line-height:1.6;color:#222;background:#fff}",
      "h1{font-size:1.6rem;border-bottom:3px solid #2C6E9B;padding-bottom:.4rem}",
      "h2{font-size:1.12rem;margin-top:2.2rem;color:#1B4B70}",
      "table{border-collapse:collapse;width:100%;margin:.7rem 0 0;font-size:.86rem}",
      "th,td{border:1px solid #d6dee5;padding:6px 9px;text-align:left;vertical-align:top}",
      "th{background:#eef4f9;font-weight:600}tr:nth-child(even) td{background:#fafcfe}",
      ".cap{font-size:.86rem;color:#555;background:#f7f9fb;border-left:4px solid #2C6E9B;",
      "padding:8px 12px;margin-bottom:.3rem}blockquote{color:#555;font-size:.9rem}",
      "</style></head><body>", "<h1>Path A — 文章级表格</h1>",
      "<p>营养-炎症分子特征 × 非小细胞肺癌免疫微环境（生成日期 2026-10-16）</p>"]

for fn, title, cap in TABLES:
    p = os.path.join(TDIR, fn)
    if not os.path.exists(p):
        continue
    rows = read_csv(p)
    md += ["## " + title, "", cap, "", md_table(rows), ""]
    ht += ["<h2>%s</h2>" % html.escape(title),
           "<div class='cap'>%s</div>" % html.escape(cap),
           html_table(rows)]

md += ["---", "",
       "## 表格使用说明", "",
       "- 主表：Table 1（基线）、Table 2（轴 × 免疫特征）、Table 3（轴 × 检查点/效应基因）、"
       "Table 4（内部验证，阴性）、Table 5（营养轴分组差异，独立）。",
       "- 附表：S1 终点敏感性、S2 特征方案敏感性、S3 综合分分组（含循环论证警示）。",
       "- **重点叙事**：Table 5 是本方案唯一不含循环论证的关联结果——营养轴与全部免疫特征零基因重叠，"
       "却与髓系/炎症通路活性系统性正相关；Table 4 如实报告预后模型无判别力。"]

ht += ["<hr><h2>表格使用说明</h2><ul>",
       "<li>主表：Table 1（基线）、Table 2（轴 × 免疫特征）、Table 3（轴 × 检查点/效应基因）、"
       "Table 4（内部验证，阴性）、Table 5（营养轴分组差异，独立）。</li>",
       "<li>附表：S1 终点敏感性、S2 特征方案敏感性、S3 综合分分组（含循环论证警示）。</li>",
       "<li><b>重点叙事</b>：Table 5 是本方案唯一不含循环论证的关联结果；"
       "Table 4 如实报告预后模型无判别力。</li></ul></body></html>"]

open(os.path.join(TDIR, "PathA_文章级表格.md"), "w", encoding="utf-8").write("\n".join(md))
open(os.path.join(TDIR, "PathA_文章级表格.html"), "w", encoding="utf-8").write("\n".join(ht))
print("OK ->", os.path.join(TDIR, "PathA_文章级表格.md"), "/ .html")
print("表格数:", sum(1 for f, _, _ in TABLES if os.path.exists(os.path.join(TDIR, f))))

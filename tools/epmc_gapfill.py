# -*- coding: utf-8 -*-
"""tools/epmc_gapfill.py — A5 定向补检（Europe PMC，中线通道）

背景：A4 的 159 条池子是按 ABSTRACT 关键词抓的，其中 T2/T3 噪声大、
T6 的检索式结构性失效（用 ABSTRACT:"external validation" 作 AND 条件，
命中的是"提到外部验证"的签名论文，而非批评它们的方法学文献）。
本脚本为 A5 的 Introduction / Discussion 做定向补检，输出可核验的候选。

用法：
  python tools/epmc_gapfill.py count          # 只报各查询命中数
  python tools/epmc_gapfill.py search         # 导出 results/refs_gapfill/*.csv

铁律：只写入 Europe PMC 实际返回的记录，不生成、不推测任何字段。
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
OUT = r"D:/WorkBuddy_Y02/out_ali_directions/pathA/results/refs_gapfill"

# 每条：key / 用途说明 / Europe PMC 检索式
QUERIES = [
    ("G1_signature_reproducibility",
     "预后基因签名不可重复 / 过拟合的方法学文献（Discussion 核心论点）",
     '(TITLE:"gene signature" OR TITLE:"gene expression signature" OR TITLE:"prognostic signature") '
     'AND (TITLE:reproducibility OR TITLE:reproducible OR TITLE:"not reproducible" '
     'OR TITLE:overfitting OR TITLE:"external validation" OR TITLE:"sample size" '
     'OR TITLE:"random validation" OR TITLE:instability)'),

    ("G2_landmark_reproducibility",
     "已知文献的定向核验（乳腺芯片年代的多随机验证与签名不唯一性）",
     '(TITLE:"multiple random validation" OR TITLE:"unique set" OR TITLE:"random gene expression signatures") '
     'AND (ABSTRACT:cancer OR ABSTRACT:outcome)'),

    ("G3_tumour_purity_confounding",
     "bulk 肿瘤转录组的纯度/间质混入混杂（限制章节的依据）",
     '(TITLE:"tumor purity" OR TITLE:"tumour purity" OR TITLE:"stromal contamination" '
     'OR TITLE:"tumor cellularity" OR TITLE:"tumour cellularity") '
     'AND (TITLE:transcriptome OR TITLE:expression OR TITLE:"gene expression" OR TITLE:microenvironment)'),

    ("G4_albumin_inflammation",
     "血清白蛋白作为负性急性期反应蛋白 / 炎症标志（营养轴的解释）",
     '(TITLE:albumin OR TITLE:hypoalbuminemia OR TITLE:hypoalbuminaemia) '
     'AND (TITLE:inflammation OR TITLE:"acute-phase" OR TITLE:"acute phase" OR TITLE:"C-reactive")'),

    ("G5_nutrition_tumour_immunity",
     "营养不良/营养指数与肿瘤免疫（营养轴 × 免疫关联的文献锚点）",
     '(TITLE:malnutrition OR TITLE:"nutritional status" OR TITLE:"nutritional index" '
     'OR TITLE:"nutritional assessment" OR TITLE:"prognostic nutritional index") '
     'AND (TITLE:immun* OR TITLE:"immune microenvironment" OR TITLE:"tumor-infiltrating" '
     'OR TITLE:"tumour-infiltrating") AND (TITLE:cancer OR TITLE:tumor OR TITLE:tumour)'),

    ("G6_hepatokine_myeloid",
     "肝源性分泌蛋白与髓系细胞活化（营养轴与巨噬/单核活化关联的机制解释）",
     '(TITLE:hepatokine* OR TITLE:"liver-derived" OR TITLE:"liver secreted" '
     'OR TITLE:"hepatic secretion") AND (ABSTRACT:macrophage OR ABSTRACT:monocyte '
     'OR ABSTRACT:myeloid OR ABSTRACT:inflammation)'),

    ("G7_tcell_infiltration_NSCLC",
     "NSCLC 中 T 细胞浸润与预后（淋巴轴方向的期望证据）",
     '(TITLE:"CD8" OR TITLE:"T cell" OR TITLE:"T-cell" OR TITLE:"tumor-infiltrating lymphocytes") '
     'AND (TITLE:"non-small cell lung" OR TITLE:NSCLC OR TITLE:"lung adenocarcinoma") '
     'AND (TITLE:prognosis OR TITLE:prognostic OR TITLE:survival)'),

    ("G8_ssgsea_vs_deconvolution",
     "通路级活性代理 vs 去卷积估计的方法学依据（2.5 节表述红线的文献支撑）",
     '(TITLE:deconvolution OR TITLE:"cell-type deconvolution" OR TITLE:"immune infiltration estimation") '
     'AND (ABSTRACT:bulk OR ABSTRACT:transcriptome OR ABSTRACT:"gene expression") '
     'AND (ABSTRACT:tumor OR ABSTRACT:tumour OR ABSTRACT:cancer)'),

    ("G9_index_to_intratumoural_immune",
     "**最接近的先行研究**：系统性营养/炎症指数与瘤内免疫特征直接关联",
     '(ABSTRACT:"prognostic nutritional index" OR ABSTRACT:"nutritional status" '
     'OR ABSTRACT:albumin OR ABSTRACT:"systemic inflammation") '
     'AND (ABSTRACT:"tumor-infiltrating" OR ABSTRACT:"tumour-infiltrating" '
     'OR ABSTRACT:"immune microenvironment" OR ABSTRACT:"tumor microenvironment") '
     'AND (ABSTRACT:NSCLC OR ABSTRACT:"lung cancer" OR ABSTRACT:"lung adenocarcinoma")'),

    ("G10_inflammation_paradox",
     "炎症的双向性（促瘤 vs 抗瘤），支撑「炎症表型决定命运」的框架",
     '(TITLE:inflammation OR TITLE:inflammatory) AND (TITLE:paradox OR TITLE:dual '
     'OR TITLE:antagonistic OR TITLE:opposing OR TITLE:heterogeneity) '
     'AND (ABSTRACT:tumor OR ABSTRACT:tumour) AND (ABSTRACT:immunity OR ABSTRACT:immunosurveillance '
     'OR ABSTRACT:"immune microenvironment")'),
]


def fetch(query, page_size):
    params = {
        "query": query,
        "format": "json",
        "pageSize": str(page_size),
        "resultType": "lite",
        "sort": "CITED desc",
    }
    url = BASE + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "pathA-REF-verify/1.0"})
    os.environ["NO_PROXY"] = "*"
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode("utf-8"))


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "count"
    if mode == "count":
        print(f"{'key':<34}{'hits':>10}  说明")
        print("-" * 110)
        for key, note, q in QUERIES:
            try:
                d = fetch(q, 1)
                print(f"{key:<34}{d['hitCount']:>10}  {note}")
            except Exception as e:
                print(f"{key:<34}{'ERR':>10}  {e}")
            time.sleep(0.4)
        return

    os.makedirs(OUT, exist_ok=True)
    import csv
    total = 0
    for key, note, q in QUERIES:
        d = fetch(q, 40)
        rows = d["resultList"]["result"]
        fn = os.path.join(OUT, f"{key}.csv")
        with open(fn, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.writer(f)
            w.writerow(["pmid", "doi", "title", "journal", "year", "pubtype",
                        "oa", "pmcid", "cited", "note"])
            for r in rows:
                w.writerow([r.get("pmid", ""), r.get("doi", ""), r.get("title", ""),
                            r.get("journalTitle", ""), r.get("pubYear", ""),
                            r.get("pubType", ""), r.get("isOpenAccess", ""),
                            r.get("pmcid", ""), r.get("citedByCount", ""), note])
        total += len(rows)
        print(f"{key:<34} hit={d['hitCount']:>8}  取回 {len(rows):>3}  -> {os.path.basename(fn)}")
        time.sleep(0.4)
    print(f"\n合计取回 {total} 条 -> {OUT}")


if __name__ == "__main__":
    main()

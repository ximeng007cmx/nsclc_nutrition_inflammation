# -*- coding: utf-8 -*-
"""A5 引用短名单：从 A4 池 + gap-fill 结果中抽取候选 PMID 的元数据，产出待核验 CSV。
只读，不改动任何源文件。输出 manuscript/refs/a5_shortlist.csv + 控制台摘要。
"""
import csv, os, glob, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOURCES = [os.path.join(ROOT, "manuscript", "refs", "03_refs.csv")] + \
          sorted(glob.glob(os.path.join(ROOT, "results", "refs_gapfill", "G*_*.csv")))
OUT = os.path.join(ROOT, "manuscript", "refs", "a5_shortlist.csv")

# 短名单（来自 gap-fill 定向检索的人工判读结果）
CAND = {
    # —— 签名/标志性研究不可重复性（T6 的真正替身）——
    "15705458": "signature reproducibility / multiple random validation",
    "15308542": "is there a unique set of genes? / reproducibility",
    "22028643": "most random signatures associate with breast cancer outcome",
    "25692416": "gene-signature reproducibility (optional)",
    # —— 最接近的前人工作（索引→瘤内免疫）——
    "25238252": "systemic inflammation + nutritional status + TIME in resected NSCLC",
    # —— 肝源性分泌蛋白 / 细胞因子 → 髓系 ——
    "33846641": "colchicine-induced hepatokines inhibit myeloid activation (Nat Metab)",
    "33436955": "LECT2 hepatokine links steatosis to inflammation via macrophages (Sci Rep)",
    # —— 肿瘤纯度混杂 ——
    "24113773": "tumour purity estimation (Aran 2013 Nat Commun)",
    "29141981": "impact of tumour purity on immune gene expression",
    # —— 瘤内 T 细胞浸润的预后意义（背景锚点）——
    "16421594": "T-cell infiltration prognostic impact",
    "26845192": "TIL / immune contexture prognosis",
    # —— 白蛋白作为负性急性期蛋白 ——
    "15660573": "albumin as negative acute-phase protein",
    "34181998": "albumin / inflammation review",
    "3584463": "albumin and inflammation (classic)",
    # —— 方法学复现硬需求（已核验）——
    "29625055": "TCGA-CDR pan-cancer clinical resource (Liu 2018 Cell)",
    "32444850": "Xena / UCSC Toil recompute (Vivian 2017 Nat Biotechnol)",
    "32055021": "Toil compendium (Vivian 2017) — 与 32444850 二者需辨明",
    "16199517": "GSVA (Hanzelmann 2013 BMC Bioinformatics)",
    "19847166": "ssGSEA (Barbie 2009 Nature)",
    "20808728": "glmnet (Friedman 2010 JSS)",
    "22960745": "TCGA-LUAD (Collisson 2014 Nature)",
    "25079552": "TCGA-LUSC (Hammerman 2012 Nature)",
    "23323831": "C-index / concordance (Harrell)",
    "25594174": "Cox model / concordance",
    "26771021": "Mayakonda / maftools? verify",
    "26209429": "R survival? verify",
    "28398314": "verify",
    "28650338": "verify",
    "29628290": "verify",
    "8668867": "verify",
}

def main():
    rows = {}
    for src in SOURCES:
        if not os.path.isfile(src):
            continue
        with open(src, encoding="utf-8-sig", newline="") as f:
            for r in csv.DictReader(f):
                pmid = (r.get("pmid") or "").strip()
                if not pmid:
                    continue
                if pmid in rows:
                    continue
                rows[pmid] = {
                    "pmid": pmid,
                    "doi": (r.get("doi") or "").strip(),
                    "title": (r.get("title") or "").strip(),
                    "journal": (r.get("journal") or "").strip(),
                    "year": (r.get("year") or "").strip(),
                    "oa": (r.get("oa") or "").strip(),
                    "cited": (r.get("cited") or "").strip(),
                    "role": CAND.get(pmid, ""),
                    "source": os.path.basename(src),
                }

    missing = [p for p in CAND if p not in rows]
    hit = [p for p in CAND if p in rows]
    out_rows = [rows[p] for p in CAND if p in rows]

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pmid", "doi", "title", "journal",
                                          "year", "oa", "cited", "role", "source"])
        w.writeheader()
        w.writerows(out_rows)

    print("候选 %d 条 → 命中 %d 条，缺 %d 条" % (len(CAND), len(hit), len(missing)))
    if missing:
        print("未在本地结果中找到的 PMID：", ", ".join(missing))
    print()
    for r in out_rows:
        print("%-9s %-5s %-28s %s" % (r["pmid"], r["year"][:4],
                                      (r["journal"] or "")[:28], (r["title"] or "")[:78]))
    print("\n→ %s" % OUT)

if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""A5 引用编号：把草稿里的 {{PMID}} 键按**首次出现顺序**替换为 [n]，并生成 Vancouver 参考文献表。
同时统计 Abstract / Introduction / Discussion 三节字数。
输入：manuscript/04_intro_discussion.draft.md + manuscript/refs/a5_refmeta.csv
输出：manuscript/04_intro_discussion.md、manuscript/refs/a5_numbering.csv
铁律：键必须全部在台账中命中；未命中的键直接报错退出，绝不静默跳过。
"""
import csv, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "manuscript", "tools"))
try:
    from qa_manuscript import english_words as _official_words
except Exception:                                    # 兜底：与 qa_manuscript 同式
    def _official_words(seg):
        seg = re.sub(r"\[[0-9,\u2013\u2014\-\s]+\]", " ", seg)
        return len(re.findall(r"[A-Za-z][A-Za-z'\-]*", seg))

DRAFT = os.path.join(ROOT, "manuscript", "04_intro_discussion.draft.md")
OUTMD = os.path.join(ROOT, "manuscript", "04_intro_discussion.md")
LEDGER = os.path.join(ROOT, "manuscript", "refs", "a5_refmeta.csv")
NUMCSV = os.path.join(ROOT, "manuscript", "refs", "a5_numbering.csv")

# NLM 期刊缩写（Europe PMC 返回全称，Vancouver 需缩写）
NLM = {
    "BMC cancer": "BMC Cancer",
    "Translational lung cancer research": "Transl Lung Cancer Res",
    "Journal of cachexia, sarcopenia and muscle": "J Cachexia Sarcopenia Muscle",
    "World journal of surgical oncology": "World J Surg Oncol",
    "Biomarkers : biochemical indicators of exposure, response, and susceptibility to chemicals": "Biomarkers",
    "Interactive cardiovascular and thoracic surgery": "Interact Cardiovasc Thorac Surg",
    "International journal of molecular sciences": "Int J Mol Sci",
    "International journal of cancer": "Int J Cancer",
    "ESMO open": "ESMO Open",
    "Cancers": "Cancers (Basel)",
    "Seminars in dialysis": "Semin Dial",
    "International journal of biological macromolecules": "Int J Biol Macromol",
    "The Journal of clinical investigation": "J Clin Invest",
    "Frontiers in endocrinology": "Front Endocrinol (Lausanne)",
    "Nature metabolism": "Nat Metab",
    "Scientific reports": "Sci Rep",
    "Hepatology (Baltimore, Md.)": "Hepatology",
    "Head & neck": "Head Neck",
    "PloS one": "PLoS One",
    "Lancet (London, England)": "Lancet",
    "Bioinformatics (Oxford, England)": "Bioinformatics",
    "PLoS computational biology": "PLoS Comput Biol",
    "Statistics in medicine": "Stat Med",
    "Nature communications": "Nat Commun",
    "Cancer immunology research": "Cancer Immunol Res",
    "British journal of cancer": "Br J Cancer",
    "Journal of thoracic oncology : official publication of the International Association for the Study of Lung Cancer": "J Thorac Oncol",
    "Immunity": "Immunity",
    "Journal of translational medicine": "J Transl Med",
    "Lung cancer (Amsterdam, Netherlands)": "Lung Cancer",
    "Computers in biology and medicine": "Comput Biol Med",
    "BMC pulmonary medicine": "BMC Pulm Med",
    "Nature": "Nature",
    "Cell": "Cell",
    "Nature biotechnology": "Nat Biotechnol",
    "Proceedings of the National Academy of Sciences of the United States of America": "Proc Natl Acad Sci U S A",
    "BMC bioinformatics": "BMC Bioinformatics",
    "Cell systems": "Cell Syst",
    "Journal of statistical software": "J Stat Softw",
    "J R Stat Soc Series B Stat Methodol": "J R Stat Soc Series B Stat Methodol",
}


def authors_vancouver(s):
    s = (s or "").strip()
    if not s:
        return ""
    parts = [p.strip() for p in s.split(",") if p.strip()]
    parts = [p.rstrip(".") for p in parts]
    if len(parts) > 6:
        return ", ".join(parts[:6]) + ", et al."
    return ", ".join(parts) + "."


def vancouver(r):
    au = authors_vancouver(r["authors"])
    j = NLM.get(r["journal"], r["journal"])
    vol = r["vol"]
    iss = r["iss"]
    pg = r["pages"]
    loc = ""
    if vol:
        loc = vol
        if iss:
            loc += "(%s)" % iss
        if pg:
            loc += ":" + pg
    elif pg:
        loc = pg
    # 标题已以 ? / ! 结尾时不再补句点（避免 "set?." 这类重复标点）
    tt = r["title"].strip().rstrip(".")
    if tt.endswith("?") or tt.endswith("!"):
        tt = tt + " "
    else:
        tt = tt + ". "
    s = "%s %s%s. %s;%s." % (au, tt, j, r["year"], loc) if loc else \
        "%s %s%s. %s." % (au, tt, j, r["year"])
    if r["doi"]:
        s += " doi:%s." % r["doi"]
    if r["pmid"] and r["pmid"] != "(none)":
        s += " PMID: %s." % r["pmid"]
    return s


def wc(text):
    return len(re.findall(r"[A-Za-z0-9][A-Za-z0-9''/%.\-–]*", text))


def main():
    ledger = {r["pmid"]: r for r in csv.DictReader(open(LEDGER, encoding="utf-8-sig"))}
    for r in csv.DictReader(open(LEDGER, encoding="utf-8-sig")):
        if r["pmid"] == "(none)":
            ledger["BENJAMINI1995"] = r
    txt = open(DRAFT, encoding="utf-8").read()

    keys = []
    for m in re.finditer(r"\{\{([A-Z0-9]+)\}\}", txt):
        k = m.group(1)
        if k not in keys:
            keys.append(k)
    missing = [k for k in keys if k not in ledger]
    if missing:
        print("❌ 草稿中以下键在台账里找不到，已中止：", ", ".join(missing))
        return 2

    num = {k: i + 1 for i, k in enumerate(keys)}
    out = re.sub(r"\{\{([A-Z0-9]+)\}\}", lambda m: "[%d]" % num[m.group(1)], txt)

    # —— 参考文献表 ——
    lines = ["", "---", "", "## References（A5 编号，按首次出现顺序；A6 与此前 13 条方法学文献合并重排）", ""]
    for k in keys:
        r = ledger[k]
        lines.append("%d. %s" % (num[k], vancouver(r)))
    lines += ["", "> 说明：编号 %d 条，全部经 Europe PMC（按 PMID 实检）与 CrossRef（按 DOI 实检）**双源核验**，"
                  "核验记录见 `refs/a5_refmeta.csv` / `a5_refmeta.json`。"
                  "Benjamini–Hochberg 1995 为 PubMed 未收录条目，仅经 CrossRef 核验。" % len(keys), ""]
    out = out.rstrip() + "\n" + "\n".join(lines)

    # —— 可选：手工维护的附录（交付核验 / 引用-内容核对表）——
    APPX = os.path.join(ROOT, "manuscript", "04_intro_discussion.appendix.md")
    if os.path.isfile(APPX):
        out = out.rstrip() + "\n\n---\n\n" + open(APPX, encoding="utf-8").read().strip() + "\n"
    with open(OUTMD, "w", encoding="utf-8") as f:
        f.write(out)

    with open(NUMCSV, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["no", "key", "pmid", "doi", "journal", "year", "group", "role"])
        for k in keys:
            r = ledger[k]
            w.writerow([num[k], k, r["pmid"], r["doi"], r["journal"], r["year"], r["group"], r["role"]])

    # —— 分节字数 ——
    def slice_between(a, b):
        i = out.find(a)
        j = out.find(b) if b else len(out)
        return out[i + len(a):j] if i >= 0 else ""

    ab = slice_between("## Abstract", "## Introduction")
    it = slice_between("## Introduction", "## Discussion")
    di = slice_between("## Discussion", "## References")
    print("引用编号 %d 条（草稿键 %d 个，全部命中）" % (len(keys), len(set(keys))))
    print("字数（**官方口径**：与 manuscript/tools/qa_manuscript.py 的 english_words() 一致，只统计以字母开头的词）：")
    print("  Abstract     %4d 词（**官方上限 <250**，Discover Oncology 指南原文 \"an abstract of less than 250 words\"）" % _official_words(ab))
    print("  Introduction %4d 词（样本实测 273–592，中位 490；自设 600–700 高于实证上限）" % _official_words(it))
    print("  Discussion   %4d 词（样本实测中位 1,376）" % _official_words(di))
    print("  两者合计     %4d 词" % (_official_words(it) + _official_words(di)))
    unused = [k for k in ledger if k not in keys]
    print()
    if unused:
        print("台账中未被 A5 正文引用的 %d 条（多为 Methods 引文，A6 并入附录 A）：" % len(unused))
        for k in unused:
            print("   %-11s %s" % (k, ledger[k]["role"][:70]))
    print("\n→ %s\n→ %s" % (OUTMD, NUMCSV))
    return 0


if __name__ == "__main__":
    sys.exit(main())

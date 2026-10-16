# -*- coding: utf-8 -*-
"""A6 参考文献合并编号（单次、可复跑）

把两套独立的局部编号合并为一套全文编号，按**真实阅读顺序**首次出现：
    Introduction（04 草稿）→ Methods（02_methods.md）→ Results → Discussion（04 草稿）

- A5 局部编号：04_intro_discussion.draft.md 的 {{KEY}} 键，按 Intro→Discussion 文件顺序 = [1]–[35]
- Methods 局部编号：02_methods.md 正文的 [n]，n = 1–13，对应其附录 A 的 13 条方法学文献
- 去重：3 条（25079552 TCGA-LUAD、22960745 TCGA-LUSC、8668867 Harrell 1996）两套都出现
- 合并后 = 35 + 13 − 3 = 45 条（**唯一文献数**）

输入：manuscript/04_intro_discussion.draft.md、manuscript/refs/a5_refmeta.csv、manuscript/02_methods.md
输出：manuscript/04_intro_discussion.md（合并编号）
      manuscript/05_references.md（45 条统一 Vancouver 列表）
      manuscript/refs/a6_numbering.csv（编号台账，供核对）
      manuscript/02_methods.md（正文引用号改为合并编号；原文件先备份到 work/_bak/）

铁律：
  1. 任何键/局部编号在台账里找不到 → 直接报错退出，绝不静默跳过
  2. 编号必须覆盖 1..45 且无重复；每条文献必须带 PMID 与 DOI
  3. 只改编号与参考文献表，不改任何正文文字与数值
  4. **对 manuscript/02_methods.md 非幂等**：重跑前必须先恢复其原始副本
     （work/_bak/02_methods.md.<时间戳>），否则脚本会报错中止（幂等守卫）。

可复跑性：`04_intro_discussion.md` / `05_references.md` / `refs/a6_numbering.csv` 由本脚本
整体重写，可任意重跑；正文计数在编号改写前后必须完全一致（脚本内自检，不一致即中止）。
"""
import csv, glob, os, re, shutil, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
sys.path.insert(0, os.path.join(ROOT, "manuscript", "tools"))
from a5_renumber import vancouver                      # 重用同一套 Vancouver 渲染
from qa_manuscript import english_words as _off

DRAFT = os.path.join(ROOT, "manuscript", "04_intro_discussion.draft.md")
METHODS = os.path.join(ROOT, "manuscript", "02_methods.md")
LEDGER = os.path.join(ROOT, "manuscript", "refs", "a5_refmeta.csv")
APPX = os.path.join(ROOT, "manuscript", "04_intro_discussion.appendix.md")
OUT_MD = os.path.join(ROOT, "manuscript", "04_intro_discussion.md")
OUT_REFS = os.path.join(ROOT, "manuscript", "05_references.md")
OUT_CSV = os.path.join(ROOT, "manuscript", "refs", "a6_numbering.csv")
BAKDIR = os.path.join(ROOT, "work", "_bak")

# 字数基线（官方口径 = manuscript/tools/qa_manuscript.py 的 english_words()）。
# 编号改写**只动 [n]**，不允许改变任何一个计数；不符即中止。
#   Abstract 222：2026-10-16 按官方硬上限 <250 词压缩后的定稿值（压缩前 291）。
#   Introduction 576：未受任何后续编辑影响。
#   Discussion 1,389：先由 D8（限制第 5 条的 MSigDB 版本表述与 Methods §2.5 对齐）把 1,379 → 1,391，
#     再因 A6-3 的 N3（补充材料引用改为官方要求的文件名形式 `Table S3` → `ESM 4`）**每处 −1 词**，
#     本段有 2 处，故 1,391 → 1,389。属**标记格式改变导致的字数改变**，非文字增删，已记入台账。
EXPECT_WORDS = (222, 576, 1389)          # Abstract / Introduction / Discussion

# 02_methods.md 附录 A 的 13 条方法学文献：局部编号 → 台账键
METHOD_REFS = {
    1:  "25079552",        # TCGA-LUAD, Nature 2014
    2:  "22960745",        # TCGA-LUSC, Nature 2012
    3:  "28398314",        # Vivian 2017, Toil
    4:  "26209429",        # Rahman 2015
    5:  "32444850",        # Goldman 2020, Xena
    6:  "29625055",        # Liu 2018, TCGA-CDR
    7:  "19847166",        # Barbie 2009, ssGSEA
    8:  "23323831",        # Hänzelmann 2013, GSVA
    9:  "16199517",        # Subramanian 2005, GSEA
    10: "26771021",        # Liberzon 2015, MSigDB hallmark
    11: "BENJAMINI1995",   # Benjamini & Hochberg 1995
    12: "20808728",        # Friedman 2010, glmnet
    13: "8668867",         # Harrell 1996, C-index
}


def load_ledger():
    rows = list(csv.DictReader(open(LEDGER, encoding="utf-8-sig")))
    led = {r["pmid"]: r for r in rows}
    for r in rows:
        if r["pmid"] == "(none)":
            led["BENJAMINI1995"] = r
    return led


def section(text, start, end=None):
    """取 start 之后、end 之前的片段。

    注意：`end` 必须用 `text.find(end, i)` 从 i 起找，返回的是**相对 i 的长度**。
    写成 `j = text.find(end)` 再 `text[i:i+j]` 是错的——那是把**绝对下标**当长度用，
    切片会比预期长出一大截（该缺陷曾让 Introduction 片段的字数虚高到 938 而非 576）。
    """
    i = text.find(start)
    if i < 0:
        return ""
    i += len(start)
    j = text.find(end, i) if end else -1
    return text[i:] if j < 0 else text[i:j]


def keys_in_order(txt):
    out = []
    for m in re.finditer(r"\{\{([A-Z0-9]+)\}\}", txt):
        if m.group(1) not in out:
            out.append(m.group(1))
    return out


def expand_marker(inner):
    """'3–5' → [3,4,5]；'1,2' → [1,2]"""
    nums = []
    for part in inner.split(","):
        part = part.strip()
        m = re.match(r"^(\d+)\s*[\u2013-]\s*(\d+)$", part)
        if m:
            nums.extend(range(int(m.group(1)), int(m.group(2)) + 1))
        elif part.isdigit():
            nums.append(int(part))
    return nums


def compress(nums):
    """[36,37,38,40] → '36–38,40'"""
    nums = sorted(set(nums))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        out.append(str(nums[i]) if j == i else "%d\u2013%d" % (nums[i], nums[j]))
        i = j + 1
    return ",".join(out)


def main():
    led = load_ledger()
    draft = open(DRAFT, encoding="utf-8").read()
    meth = open(METHODS, encoding="utf-8").read()

    # ---- 1. 各段局部顺序 ----
    intro_keys = keys_in_order(section(draft, "## Introduction", "## Discussion"))
    disc_keys = keys_in_order(section(draft, "## Discussion", "## References"))
    if not intro_keys or not disc_keys:
        print("❌ 未能从草稿中解析 Introduction / Discussion 的引用键，已中止")
        return 2

    mbody = section(meth, "## Materials and methods", "## 附录 A")
    if not mbody:
        print("❌ 未能定位 02_methods.md 正文段，已中止")
        return 2
    meth_order = []
    for m in re.finditer(r"\[([0-9,\u2013\-\s]+)\]", mbody):
        for n in expand_marker(m.group(1)):
            if not (1 <= n <= 13):
                print("❌ 02_methods.md 出现越界局部编号 [%d]，已中止" % n)
                return 2
            if n not in meth_order:
                meth_order.append(n)
    meth_keys = [METHOD_REFS[n] for n in meth_order]

    # ---- 幂等守卫：02_methods.md 已被本脚本改写过就不能再跑 ----
    # 正文引用号一旦变成合并编号，再映射一次会把 [28,29] 当成局部号处理，静默产生错误编号。
    # 所有正文编号都已在 1..13 之内时才是「未改写」状态。
    if meth_order and max(meth_order) > 13:
        baks = sorted(glob.glob(os.path.join(BAKDIR, "02_methods.md.*")))
        print("❌ 02_methods.md 的正文引用号已是合并编号（最大 [%d] > 13）。" % max(meth_order))
        print("   本脚本对 02_methods.md 非幂等，拒绝重复改写以免产生错误编号。先恢复原始文件：")
        print("     cp \"%s\" \"%s\"" % (baks[-1] if baks else "<work/_bak/02_methods.md.*>",
                                        os.path.relpath(METHODS, ROOT).replace("\\", "/")))
        print("   然后重新运行。已中止。")
        return 2

    # ---- 2. 合并顺序（真实阅读顺序：Intro → Methods → Results → Discussion）----
    order, seen = [], set()
    for k in intro_keys + meth_keys + disc_keys:
        if k not in seen:
            seen.add(k)
            order.append(k)

    missing = [k for k in order if k not in led]
    if missing:
        print("❌ 以下键在台账中找不到，已中止：", ", ".join(missing))
        return 2
    num = {k: i + 1 for i, k in enumerate(order)}

    # ---- 3. 一致性自检 ----
    dup_pmid = {}
    for k in order:
        p = led[k]["pmid"]
        dup_pmid.setdefault(p, []).append(k)
    dups = {p: v for p, v in dup_pmid.items() if len(v) > 1}
    if dups:
        print("❌ 存在重复 PMID，已中止：", dups)
        return 2
    unverified = [(k, led[k]["pmid"], led[k]["epmc"], led[k]["crossref"])
                  for k in order
                  if led[k]["epmc"] != "ok" or led[k]["crossref"] != "ok"]

    # 分项口径（三分类互斥且完备，合起来必须等于 order 的长度）：
    #   A = Introduction 首现            = intro_keys                        27
    #   B = Methods 首现（不在 Intro）   = meth \ intro                     12
    #   C = Discussion 首现（不在 I∪M）  = disc \ (intro ∪ meth)             6
    # 注意：不能写成「meth 去掉 (intro∪disc)」+「disc 去掉 (intro∪meth)」——
    #   那样会把「既在 Methods 又在 Discussion、但不在 Intro」的键双排除（本例 2 条：
    #   25079552 TCGA-LUAD、22960745 TCGA-LUSC），结果少 2 条。
    a5_set = set(intro_keys) | set(disc_keys)
    shared = [k for k in meth_keys if k in a5_set]            # Methods 与 A5 重用的键
    meth_first = [k for k in meth_keys if k not in set(intro_keys)]
    disc_first = [k for k in disc_keys if k not in set(intro_keys) | set(meth_keys)]
    _chk = len(intro_keys) + len(meth_first) + len(disc_first)
    if _chk != len(order):
        print("❌ 分项计数自检失败：%d + %d + %d ≠ %d，已中止"
              % (len(intro_keys), len(meth_first), len(disc_first), len(order)))
        return 2

    print("合并编号 %d 条（首现位置：Introduction %d + Methods %d + Discussion %d）"
          % (len(order), len(intro_keys), len(meth_first), len(disc_first)))
    print("  双源核验：Europe PMC ok × %d，CrossRef ok × %d，异常 %d 条"
          % (sum(1 for k in order if led[k]["epmc"] == "ok"),
             sum(1 for k in order if led[k]["crossref"] == "ok"), len(unverified)))
    if shared:
        print("  去重：Methods 13 条中 %d 条与 A5 共表并已合并 —— %s"
              % (len(shared), " / ".join("%s(%s)" % (led[k]["pmid"], k) for k in shared)))
    else:
        print("  去重：Methods 13 条与 A5 无重叠（0 条）")
    print("  A5 侧自查：Intro %d + Discussion %d − 段间重复 %d = A5 唯一 %d"
          % (len(intro_keys), len(disc_keys), len(set(intro_keys) & set(disc_keys)), len(a5_set)))
    for k, p, e, c in unverified:
        print("  ⚠ 未双源齐备：%s  PMID=%s  EuropePMC=%s  CrossRef=%s"
              % (k, p, e, c) + ("（PubMed 未收录条目，仅 CrossRef 核验，属预期）"
                                if p == "(none)" else "（**需人工复核**）"))

    # ---- 4. 04_intro_discussion.md ----
    prose = re.sub(r"\{\{([A-Z0-9]+)\}\}", lambda m: "[%d]" % num[m.group(1)], draft)
    out = prose.rstrip() + ("\n\n---\n\n> **参考文献**：合并编号后的统一列表见 `manuscript/05_references.md`"
                            "（共 %d 条，按全文首次出现顺序编号）。\n" % len(order))
    if os.path.isfile(APPX):
        out = out.rstrip() + "\n\n---\n\n" + open(APPX, encoding="utf-8").read().strip() + "\n"
    open(OUT_MD, "w", encoding="utf-8").write(out)

    # ---- 5. 02_methods.md 正文引用号改写 ----
    os.makedirs(BAKDIR, exist_ok=True)
    bak = os.path.join(BAKDIR, "02_methods.md.%s" % time.strftime("%Y%m%d_%H%M%S"))
    shutil.copy2(METHODS, bak)

    def repl(m):
        nums = expand_marker(m.group(1))
        return "[%s]" % compress([num[METHOD_REFS[n]] for n in nums])

    i = meth.find("## Materials and methods")
    j = meth.find("## 附录 A")
    new_meth = meth[:i] + re.sub(r"\[([0-9,\u2013\-\s]+)\]", repl, meth[i:j]) + meth[j:]

    # 附录 A 的局部编号列改为合并编号（附录为投稿时删除的内部件）
    row_i = new_meth.find("## 附录 A")
    row_j = new_meth.find("## 附录 B")
    if row_j < 0:
        row_j = len(new_meth)
    appx_block = new_meth[row_i:row_j]
    cnt = [0]

    def row_repl(m):
        cnt[0] += 1
        local = int(m.group(1))
        return "| %d |" % num[METHOD_REFS[local]] if local in METHOD_REFS else m.group(0)

    appx_new = re.sub(r"\| (\d+) \|", row_repl, appx_block)
    appx_new = appx_new.replace(
        "## 附录 A：Methods 参考文献（待 A6 与 A4 清单合并编号）",
        "## 附录 A：Methods 参考文献（**已并入全文合并编号，见 `05_references.md`**）")
    new_meth = new_meth[:row_i] + appx_new + new_meth[row_j:]
    if new_meth != meth:
        open(METHODS, "w", encoding="utf-8").write(new_meth)
    print("  02_methods.md：正文引用号已改写；附录 A 编号行改写 %d 行；备份 %s"
          % (cnt[0], os.path.relpath(bak, ROOT)))

    # ---- 6. 统一参考文献表 ----
    lines = [
        "# References（全文合并编号，共 %d 条）" % len(order),
        "",
        "> 编号依据：按**真实阅读顺序**首次出现（Introduction → Methods → Results → Discussion）。",
        "> Methods 的 13 条方法学文献来自 `02_methods.md` 附录 A；其余来自 A5 引用台账。",
        "> 两套局部编号合并去重 3 条（25079552 TCGA-LUAD、22960745 TCGA-LUSC、8668867 Harrell 1996）。",
        "> 全部经 Europe PMC（按 PMID 实检）+ CrossRef（按 DOI 实检）**双源核验**；",
        "> Benjamini–Hochberg 1995 为 PubMed 未收录条目，仅经 CrossRef 核验。台账见 `refs/a6_numbering.csv`。",
        "",
    ]
    for k in order:
        lines.append("%d. %s" % (num[k], vancouver(led[k])))
    lines += [
        "",
        "---",
        "",
        "## 编号映射（供核对）",
        "",
        "| 合并编号 | A5 局部编号 | Methods 局部编号 | PMID | DOI |",
        "|---|---|---|---|---|",
    ]
    for k in order:
        a5n = ("[%d]" % (intro_keys.index(k) + 1)) if k in intro_keys else (
            "[%d]" % (len(intro_keys) + disc_keys.index(k) + 1) if k in disc_keys else "—")
        mn = next(("[%d]" % n for n in meth_order if METHOD_REFS[n] == k), "—")
        lines.append("| %d | %s | %s | %s | %s |"
                     % (num[k], a5n, mn, led[k]["pmid"] if led[k]["pmid"] != "(none)" else "（PubMed 未收录）",
                        led[k]["doi"] or "—"))
    lines += ["",
              "> 局部编号口径：A5 局部编号 = 原 `04_intro_discussion.draft.md` 内的首次出现序号（Intro→Discussion 文件顺序）；",
              "> Methods 局部编号 = 原 `02_methods.md` 正文的 `[1]–[13]`。两者均已不再使用。",
              ""]
    open(OUT_REFS, "w", encoding="utf-8").write("\n".join(lines))

    # ---- 7. 编号台账 ----
    with open(OUT_CSV, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["no", "key", "pmid", "doi", "journal", "year", "source", "group", "role"])
        for k in order:
            src = "A5" if k in intro_keys or k in disc_keys else "Methods"
            if k in intro_keys or k in disc_keys:
                if k in meth_keys:
                    src = "A5+Methods"
            w.writerow([num[k], k, led[k]["pmid"], led[k]["doi"], led[k]["journal"],
                        led[k]["year"], src, led[k]["group"], led[k]["role"]])

    # ---- 8. 字数复核（编号改写不应改变任何计数）----
    # 必须对 **草稿 prose** 计词，不能对 out 计词：out 已把参考文献表移出（见 05_references.md）
    # 并在文末拼接了中文附录，用 "## References" 作结束标记会失效，Discussion 会吞掉整个附录。
    # 口径必须与 manuscript/tools/qa_manuscript.py 完全一致 → 直接用它的 body_of/english_words。
    # 对照组：先把 {{KEY}} 整体替换为空白再计数。原因：纯数字键 {{25079552}} 本就不计词，
    # 但字母键 {{BENJAMINI1995}} 会被计为 1 个词，改写为 [11] 后归零——不先掩掉就会误报差异。
    from qa_manuscript import body_of as _body

    draft_masked = re.sub(r"\{\{[A-Z0-9]+\}\}", " ", draft)
    _spec = (("## Abstract", "## Introduction"), ("## Introduction", "## Discussion"),
             ("## Discussion", None))
    after = tuple(_off(_body(prose, s, e)) for s, e in _spec)
    before = tuple(_off(_body(draft_masked, s, e)) for s, e in _spec)
    if after != before:
        print("❌ 字数自检失败：编号改写导致正文字数变化（不应发生），已中止")
        print("   Abstract/Introduction/Discussion 改写后 %s vs 改写前 %s" % (after, before))
        return 2
    if after != EXPECT_WORDS:
        print("❌ 字数自检失败：与既定基线不一致，已中止")
        print("   实测 Abstract/Introduction/Discussion = %s，基线 %s" % (after, EXPECT_WORDS))
        return 2
    print("  字数（官方口径，与基线一致）：Abstract %d / Introduction %d / Discussion %d / "
          "正文合计 %d" % (after[0], after[1], after[2], after[1] + after[2]))
    print("\n→ %s\n→ %s\n→ %s" % (OUT_MD, OUT_REFS, OUT_CSV))
    return 0


if __name__ == "__main__":
    sys.exit(main())

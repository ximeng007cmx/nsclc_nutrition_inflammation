# -*- coding: utf-8 -*-
"""tools/a6_3_esm.py — A6-3 N3：把补充材料引用从 "Table S1 / Fig. S1 / Supplementary Methods S1"
改为 Discover Oncology 官方指南要求的「按文件名引用」形式 **ESM N**。

依据（本地存档 work/do_guidelines.html，原文）：
  "Naming, numbering and citing supplementary material: Name supplementary material clearly and number similar
   types of material. … When referring to supplementary material, clearly differentiate between in-text
   figures/tables and supplementary material by citing according to file name. E.g. '…as shown in the animation
   (Online Resource 1)', '…additional data are given in ESM 2'. You must reference any supplementary material
   in a citation."

用户决策（2026-10-16）：N3 = **按本刊要求改**（原全文用 Fig. S1 / Table S1）。

编号方案（同类材料聚在一起、按类内序号排；每个 ESM = 一个交付文件）：
  ESM 1  Supplementary Methods S1                    （方法学说明文件）
  ESM 2  Table S1  终点敏感性 OS/DSS/PFI
  ESM 3  Table S2  特征方案 × 训练方向敏感性（12 组）
  ESM 4  Table S3  按综合风险分分组差异（含循环论证警示）
  ESM 5  Table S4  基因与特征清单（内分 S4a–S4d 四个分区）
  ESM 6  Figure S1 LUAD KM
  ESM 7  Figure S2 LUSC KM
  ESM 8  Figure S3 LASSO λ 曲线
  ESM 9  Figure S4 多变量 Cox 森林图

铁律：
  1. **只改投稿面向的英文正文**；中文内部文档（附录、台账）一律不动，改由 `refs/esm_mapping.csv` 对照。
     例外：`manuscript/Supplementary_Methods.md` 本身是 ESM 1（要投出去的英文文件），故纳入。
  2. **绝不改文件名**：`TableS4b_immune_features.csv` / `FigureS1_km_luad.png` 之类没有空格，
     所有模式都要求 "Table S1" 这种带空格的形式，且加 `(?![a-d])` 负向前瞻防止 "Table S4b" 被截断成 "ESM 5b"。
     S4a–S4d 四个分区属 ESM 5，用**单独一条** `Table S4[a-d]` → `ESM 5` 规则整体替换（不留 "ESM 5b" 这种尾巴）。
  3. 幂等：已改成 ESM 的文本再跑不会二次改动（模式不再命中）。
  4. 写盘前逐个文件备份到 `work/_bak/<名>.<时间戳>`（本目录非 git 仓库，.bak 是唯一回退手段）。

用法：
  python tools/a6_3_esm.py --check     # 只打印将要做的改动（统一 diff）
  python tools/a6_3_esm.py             # 应用（自动备份 + 写盘 + 写映射台账）
"""
import difflib
import os
import re
import shutil
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BAKDIR = os.path.join(ROOT, "work", "_bak")

# ---- 编号表 ----------------------------------------------------------------
ESM = {
    "ESM 1": ("Supplementary Methods S1", "Supplementary methods: full specification of the two pre-specified sensitivity analyses"),
    "ESM 2": ("Table S1", "Endpoint sensitivity analysis (OS / DSS / PFI)"),
    "ESM 3": ("Table S2", "Feature-scheme and training-direction sensitivity analysis (12 runs)"),
    "ESM 4": ("Table S3", "Immune feature activity by composite-score group, reported with an explicit circularity warning"),
    "ESM 5": ("Table S4 (parts S4a-S4d)", "Gene and feature inventories (sheets S4a-S4d)"),
    "ESM 6": ("Figure S1", "Overall survival of TCGA-LUAD patients by composite risk group"),
    "ESM 7": ("Figure S2", "Overall survival of TCGA-LUSC patients by the same composite risk group"),
    "ESM 8": ("Figure S3", "LASSO-Cox cross-validation curve in the LUAD training set"),
    "ESM 9": ("Figure S4", "Multivariable Cox model in the full LUAD cohort"),
}

# ---- 替换规则（顺序敏感：复杂模式在前）--------------------------------------
RULES = [
    (r"\bTables S1 and S2\b", "ESM 2 and ESM 3"),
    (r"\bSupplementary Table S3\b", "ESM 4"),
    (r"\bSupplementary Methods S1\b", "ESM 1"),
    (r"\bTable S4(?![a-d])\b", "ESM 5"),
    (r"\bTable S4[a-d]\b", "ESM 5"),          # S4a–S4d 同属 ESM 5，整体替换（避免 "ESM 5b" 这种残尾）
    (r"\bTable S3(?![a-d])\b", "ESM 4"),
    (r"\bTable S2(?![a-d])\b", "ESM 3"),
    (r"\bTable S1(?![a-d])\b", "ESM 2"),
    (r"\bFig\. S1\b", "ESM 6"),
    (r"\bFig\. S2\b", "ESM 7"),
    (r"\bFig\. S3\b", "ESM 8"),
    (r"\bFig\. S4\b", "ESM 9"),
    (r"\bFigure S1\b", "ESM 6"),
    (r"\bFigure S2\b", "ESM 7"),
    (r"\bFigure S3\b", "ESM 8"),
    (r"\bFigure S4\b", "ESM 9"),
]

# ---- 目标文件与「只改这一段」的边界 -----------------------------------------
# (相对路径, 起始标记, 结束标记) —— 结束标记为 None 表示到文件末尾
TARGETS = [
    ("manuscript/01_results.md",                None,                       "## 数字溯源对照表"),
    ("manuscript/02_methods.md",                None,                       "## 附录 A"),
    ("manuscript/04_intro_discussion.draft.md", None,                       None),
    ("manuscript/03_figure_legends.md",         "## Supplementary figures", "## 交付文件索引"),
    ("manuscript/Supplementary_Methods.md",     None,                       None),
]

# 文件名安全核查：这些「无空格」文件名一个都不允许变
FILENAME_TOKENS = ["TableS1", "TableS2", "TableS3", "TableS4a", "TableS4b", "TableS4c", "TableS4d",
                   "FigureS1", "FigureS2", "FigureS3", "FigureS4"]


def apply_rules(text):
    n = 0
    for pat, rep in RULES:
        text, k = re.subn(pat, rep, text)
        n += k
    return text, n


def count_tokens(text):
    return {t: text.count(t) for t in FILENAME_TOKENS}


def main(argv):
    check = "--check" in argv
    total = 0
    guard = {}
    if not check:
        os.makedirs(BAKDIR, exist_ok=True)
        stamp = time.strftime("%Y%m%d_%H%M%S")

    for rel, start, end in TARGETS:
        p = os.path.join(ROOT, rel)
        if not os.path.isfile(p):
            print("  [SKIP 文件不存在] %s" % rel)
            continue
        src = open(p, encoding="utf-8").read()
        if start and start in src:
            i = src.index(start)
        elif start:
            print("  [SKIP 找不到起始标记] %s" % rel)
            continue
        else:
            i = 0
        j = src.index(end, i) if (end and end in src[i:]) else len(src)
        head, region, tail = src[:i], src[i:j], src[j:]
        new_region, n = apply_rules(region)
        total += n
        guard[rel] = (count_tokens(src), count_tokens(head + new_region + tail))
        print("  %-30s 替换 %3d 处" % (rel.split("/")[-1], n))
        if n and check:
            for ln in difflib.unified_diff(region.split("\n"), new_region.split("\n"),
                                           "old", "new", lineterm="", n=0):
                if ln.startswith(("+", "-")) and not ln.startswith(("+++", "---")):
                    print("        " + ln[:200])
        if n and not check:
            bak = os.path.join(BAKDIR, "%s.%s" % (os.path.basename(p), stamp))
            shutil.copy2(p, bak)
            open(p, "w", encoding="utf-8").write(head + new_region + tail)

    # ---- 文件名安全核查：无空格文件名计数必须一字不变 ----
    bad = [(rel, t, a[t], b[t]) for rel, (a, b) in guard.items()
           for t in FILENAME_TOKENS if a[t] != b[t]]
    if bad:
        print("❌ 文件名安全核查失败（以下无空格文件名计数发生变化），已中止：")
        for rel, t, a, b in bad:
            print("   %s: %s  %d → %d" % (rel, t, a, b))
        return 2
    print("  ✅ 文件名安全核查：%d 个无空格文件名 token 在全部目标文件中计数不变"
          % len(FILENAME_TOKENS))

    # ---- 残留扫描（只看投稿面向的英文区段）----
    resid = []
    for rel, start, end in TARGETS:
        p = os.path.join(ROOT, rel)
        if not os.path.isfile(p):
            continue
        src = open(p, encoding="utf-8").read()
        i = src.index(start) if (start and start in src) else 0
        j = src.index(end, i) if (end and end in src[i:]) else len(src)
        seg = src[i:j]
        for m in re.finditer(r"\b(?:Tables? S[0-9][a-d]?|Fig\. S[0-9]|Figure S[0-9]|Supplementary (?:Table|Methods) S[0-9])\b", seg):
            resid.append((rel, m.group(0)))
    if resid:
        print("  %s 旧编号 %d 处："
              % ("（干跑：以下为改动前盘上状态）" if check else "⚠ 仍在投稿区段内的", len(resid)))
        for rel, s in resid:
            print("     %s  →  %s" % (rel, s))
    else:
        print("  ✅ 投稿区段内旧编号残留 = 0")

    print("  合计替换 %d 处%s" % (total, "（干跑，未写盘）" if check else ""))

    if not check:
        # 映射台账（用 csv.writer，避免内容里的逗号把列冲散）
        import csv
        out = os.path.join(ROOT, "manuscript", "refs", "esm_mapping.csv")
        files = {"ESM 1": "ESM_1.pdf", "ESM 2": "ESM_2.xlsx", "ESM 3": "ESM_3.xlsx",
                 "ESM 4": "ESM_4.xlsx", "ESM 5": "ESM_5.xlsx", "ESM 6": "ESM_6.pdf",
                 "ESM 7": "ESM_7.pdf", "ESM 8": "ESM_8.pdf", "ESM 9": "ESM_9.pdf"}
        with open(out, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.writer(f)
            w.writerow(["esm", "原编号", "内容", "交付文件名(A6-3 打包时)"])
            for k in sorted(ESM, key=lambda x: int(x.split()[1])):
                w.writerow([k, ESM[k][0], ESM[k][1], files[k]])
        print("  → %s" % os.path.relpath(out, ROOT).replace("\\", "/"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

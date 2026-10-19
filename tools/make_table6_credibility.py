# -*- coding: utf-8 -*-
"""
make_table6_credibility.py — 生成主文 Table 6（阴性结果可信度量表）。

数据全部实测取自 results/*.tsv（2026-10-18 R/08 重跑后口径），不手打：
  axis_univariate_cox_luad.tsv   → 行 1：lymphocyte 单轴 HR/P
  validation_summary.tsv         → 行 1：全队列 in-sample C（LUAD train）
  geneset_sensitivity.tsv        → 行 1/5：S0/S1/S2 训练 C 与 LUSC 验证 C
  internal_repeated_cv.tsv       → 行 3：100 次分割 held-out C 的 type-7 分位数、n_sel 计数
  random_panel_null_summary.tsv  → 行 4：单轴 C vs expr_matched 随机面板（p_emp）
  incremental_cox_terms.tsv      → 行 6：M1/M4 risk_score HR
  incremental_cox_tests.tsv      → 行 6：M5→M4 LR P
  external_validation.tsv        → 行 7：复合分外部 C/CI、GSE31210 neutrophil HR/FDR

设计约束：
  * 表内英文一律不用 em dash（沿用正文红线）；区间用 " to "。
  * 每个格式化数值都对照硬编码期望值断言，数值漂移即报错退出（同 check_results_numbers 家法）。
"""
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = RES / "tables" / "Table6_negative_result_credibility.csv"


def read_tsv(name):
    with open(RES / name, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def f2(x, nd=3):
    return f"{float(x):.{nd}f}"


def q7(vals, p):
    """R quantile type-7（与 independent_recompute.py 一致）。"""
    s = sorted(vals)
    h = (len(s) - 1) * p
    lo = math.floor(h)
    hi = math.ceil(h)
    return s[lo] + (h - lo) * (s[hi] - s[lo])


def expect(got, want, tag):
    if got != want:
        print(f"[Table6] 数值漂移 {tag}: got {got!r}, want {want!r}", file=sys.stderr)
        sys.exit(1)


# ---------- 行 1：模型有能力拟合信号 ----------
uni = {r["axis"]: r for r in read_tsv("axis_univariate_cox_luad.tsv")}
lym_hr = f2(uni["lymphocyte"]["HR"])
lym_p = f2(uni["lymphocyte"]["p"])
expect(lym_hr, "0.843", "row1 HR")
expect(lym_p, "0.025", "row1 P")
val = {r["cohort"]: r for r in read_tsv("validation_summary.tsv")}
c_insample = f2(val["TCGA-LUAD (train)"]["cindex"])
expect(c_insample, "0.555", "row1 in-sample C")
gs_all = read_tsv("geneset_sensitivity.tsv")
gs = {}
for r in gs_all:
    if r["endpoint"] == "OS" and r["test"] == "LUSC":
        assert r["scheme"] not in gs, "endpoint=OS 行重复"
        gs[r["scheme"]] = r
assert len(gs) == 3, f"期望 3 个 OS scheme，得 {sorted(gs)}"
c_big = f2(gs["S2_genes48"]["c_train"])
expect(c_big, "0.681", "row1 big-scheme C")

# ---------- 行 3：100 次随机分割 ----------
rcv = read_tsv("internal_repeated_cv.tsv")
ctests = [float(r["c_test"]) for r in rcv]
n_single = sum(1 for r in rcv if r["n_sel"].strip() == "1")
med = f2(q7(ctests, 0.50))
lo = f2(q7(ctests, 0.025))
hi = f2(q7(ctests, 0.975))
expect(med, "0.542", "row3 median")
expect(lo, "0.456", "row3 q025")
expect(hi, "0.600", "row3 q975")
expect(str(n_single), "92", "row3 n_single")

# ---------- 行 4：随机面板 null ----------
rnd = {}
for r in read_tsv("random_panel_null_summary.tsv"):
    if r["null_type"] == "expr_matched" and r["metric"] == "os_Cindex_luad":
        rnd.setdefault(r["axis"], r)
c_nut = f2(rnd["nutrition_liver"]["obs"])
c_neu = f2(rnd["neutrophil"]["obs"])
c_lym = f2(rnd["lymphocyte"]["obs"])
p_nut = f2(rnd["nutrition_liver"]["p_emp"], 2)
p_neu = f2(rnd["neutrophil"]["p_emp"], 2)
p_lym = f2(rnd["lymphocyte"]["p_emp"], 2)
p_lo = min(p_nut, p_neu, p_lym)
p_hi = max(p_nut, p_neu, p_lym)
c_min = min(c_nut, c_neu, c_lym)
c_max = max(c_nut, c_neu, c_lym)
expect(c_nut, "0.535", "row4 C nutrition")
expect(c_neu, "0.470", "row4 C neutrophil")
expect(c_lym, "0.445", "row4 C lymphocyte")
expect(p_nut, "0.08", "row4 P nutrition")
expect(p_neu, "0.74", "row4 P neutrophil")
expect(p_lym, "0.90", "row4 P lymphocyte")

# ---------- 行 5：敏感性（过拟合签名） ----------
# S0 的 c_train 与 validation_summary 的全队列 in-sample C 是同一量，但本 tsv 只存
# 4 位预修约字符串 "0.5545"（二进制浮点会修约成 0.554）；取全精度源与正文口径一致。
s1c, s2c = (f2(gs[s]["c_train"]) for s in ("S1_msigdb6", "S2_genes48"))
s0c = c_insample
s0t, s1t, s2t = (f2(gs[s]["c_test"]) for s in ("S0_prior16", "S1_msigdb6", "S2_genes48"))
expect(s0c, "0.555", "row5 S0 train"); expect(s1c, "0.583", "row5 S1 train")
expect(s2c, "0.681", "row5 S2 train"); expect(s0t, "0.501", "row5 S0 test")
expect(s1t, "0.518", "row5 S1 test");  expect(s2t, "0.482", "row5 S2 test")

# ---------- 行 6：估计器/模型族独立 ----------
terms = read_tsv("incremental_cox_terms.tsv")
hr_m1 = f2(next(r["HR"] for r in terms if r["model"] == "M1 clinical + risk_score" and r["term"] == "risk_score"), 2)
hr_m4 = f2(next(r["HR"] for r in terms if r["model"] == "M4 M1 + purity + ImmuneScore" and r["term"] == "risk_score"), 2)
expect(hr_m1, "15.49", "row6 HR M1")
expect(hr_m4, "2.51", "row6 HR M4")
tests = read_tsv("incremental_cox_tests.tsv")
lr_m54 = f2(next(r["p"] for r in tests if r["comparison"].startswith("M5 -> M4 (+risk_score")), 3)
expect(lr_m54, "0.738", "row6 LR P")

# ---------- 行 7：外部冻结模型 ----------
ext = read_tsv("external_validation.tsv")
comp = {r["cohort"]: r for r in ext if r["analysis"].startswith("composite")}
axes = {(r["cohort"], r["analysis"]): r
        for r in ext if r["analysis"].startswith("single axis")}
c684 = f2(comp["GSE68465"]["cindex"])
lo684 = f2(comp["GSE68465"]["lo"])
hi684 = f2(comp["GSE68465"]["hi"])
c312 = f2(comp["GSE31210"]["cindex"])
lo312 = f2(comp["GSE31210"]["lo"])
hi312 = f2(comp["GSE31210"]["hi"])
neu312 = axes[("GSE31210", "single axis: neutrophil")]
neu_hr = f2(neu312["HR"])
neu_fdr = f2(neu312["fdr"], 4)
expect(c684, "0.473", "row7 C GSE68465")
expect(lo684, "0.431", "row7 lo GSE68465")
expect(hi684, "0.517", "row7 hi GSE68465")
expect(c312, "0.593", "row7 C GSE31210")
expect(lo312, "0.487", "row7 lo GSE31210")
expect(hi312, "0.698", "row7 hi GSE31210")
expect(neu_hr, "1.516", "row7 neutrophil HR")
expect(neu_fdr, "0.0094", "row7 neutrophil FDR")

# ---------- 组装 ----------
rows = [
    [
        "Line of evidence",
        "Result",
        "Alternative explanation it rules out",
    ],
    [
        "Capability of the analysis to detect a signal",
        f"Cross-validated screening on the full LUAD cohort retained the lymphocyte axis (univariate HR {lym_hr}, P = {lym_p}); in-sample C-index {c_insample}, rising to {c_big} with larger feature schemes",
        "That the model was unable to fit any structure in these data",
    ],
    [
        "Design against information leakage",
        "Feature screening, standardisation and the risk cut-point were fixed within the training data, and evaluation used held-out samples only",
        "Optimism from fitting and evaluating performance on the same samples",
    ],
    [
        "Internal replication over 100 random 7:3 splits",
        f"Median held-out C-index {med} (2.5th to 97.5th percentile {lo} to {hi}); the predictor set was a single axis in {n_single} of 100 splits",
        "That the null result reflects one unlucky split",
    ],
    [
        "Expression-matched random-panel null",
        f"All three single-axis C-indices ({c_min} to {c_max}) lay within their random-panel distributions (empirical P = {p_lo} to {p_hi})",
        "That the composite carries prognostic information that random expression-matched panels do not",
    ],
    [
        "Sensitivity analyses",
        f"As features were added the training C-index rose ({s0c} to {s1c} to {s2c}) while the cross-histology validation C-index fell ({s0t} to {s1t} to {s2t}); no change of endpoint or cross-histology validation was significant",
        "That the null is an artefact of one modelling choice",
    ],
    [
        "Estimator and model-family independence",
        f"An independent marker-based estimator reproduced the association directions, and in multivariable models the composite attenuated from HR {hr_m1} to {hr_m4} once immune content entered (likelihood-ratio P = {lr_m54}) with no gain in discrimination",
        "That the null depends on the scoring method or the survival model family",
    ],
    [
        "External application of the frozen model",
        f"C-index {c684} (95% CI {lo684} to {hi684}) in GSE68465 and {c312} ({lo312} to {hi312}) in GSE31210; the only FDR-significant single-axis association (GSE31210 neutrophil, HR {neu_hr}, FDR = {neu_fdr}) ran opposite to TCGA",
        "That the null is specific to the TCGA compendium",
    ],
]

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w", newline="", encoding="utf-8") as f:
    csv.writer(f).writerows(rows)

# 表内红线自查：英文列不得含 em dash
bad = [c for r in rows for c in r if "\u2014" in c]
if bad:
    print(f"[Table6] 表内出现 em dash: {bad}", file=sys.stderr)
    sys.exit(1)

print(f"[Table6] OK -> {OUT} ({OUT.stat().st_size} B, {len(rows) - 1} evidence rows)")

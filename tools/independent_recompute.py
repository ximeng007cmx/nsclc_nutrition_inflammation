#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Independent recomputation of every headline number in the manuscript.

Reads ONLY the raw result tables under results/ and recomputes the statistics
from scratch (no reuse of numbers written in the manuscript or in Table 4 CSV).
Each block prints: recomputed value | manuscript claim | verdict.

Files used:
  results/table4_nutrition_axis_raw.tsv        (nutrition-axis median split)
  results/immune_axis_correlation.tsv          (39 axis-feature pairs)
  results/checkpoint_correlation.tsv           (57 axis-transcript pairs)
  results/axis_univariate_cox_luad.tsv
  results/model_coef.tsv
  results/internal_split_summary.tsv
  results/internal_repeated_cv.tsv
  results/_audit_internal_selected.tsv         (axis identity per split)
  results/multivariate_cox_luad.tsv
  results/geneset_sensitivity.tsv
  results/validation_summary.tsv
"""
import csv
import os
import statistics as st

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(BASE, "results")


def load(name, delim="\t"):
    with open(os.path.join(R, name), encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter=delim))


def f(x):
    return float(x)


checks = []


def chk(label, got, claim, tol=0.0015, fmt="%.3f"):
    """Numeric comparison with display rounding tolerance."""
    if isinstance(got, str):
        ok = (got == claim)
        print("%-52s %-28s %-24s %s" % (label, got, claim, "OK" if ok else "**DIFF**"))
        checks.append(ok)
        return
    ok = abs(got - claim) <= tol
    print("%-52s %-28s %-24s %s" % (label, fmt % got, fmt % claim, "OK" if ok else "**DIFF**"))
    checks.append(ok)


print("=" * 108)
print("A. Table 4 — nutrition-axis median split (results/table4_nutrition_axis_raw.tsv)")
print("=" * 108)
t4 = load("table4_nutrition_axis_raw.tsv")
print("rows: %d | n_high set: %s | n_low set: %s" % (
    len(t4), sorted({r["n_high"] for r in t4}), sorted({r["n_low"] for r in t4})))
diffs = []
for r in sorted(t4, key=lambda x: -f(x["diff"])):
    d = f(r["median_high"]) - f(r["median_low"])
    diffs.append(d)
    print("  %-46s high=%.3f low=%.3f diff=%+.4f  (file diff %+.4f, p=%s)"
          % (r["immune_feature"], f(r["median_high"]), f(r["median_low"]), d, f(r["diff"]), r["p"]))
print("  -> all diffs positive:", all(d > 0 for d in diffs), "| max %.3f | min %.3f" % (max(diffs), min(diffs)))

chk("Table 4 max diff (monocyte activation)", max(diffs), 0.067)
chk("Table 4 min diff (neutrophil degranulation)", min(diffs), 0.009)
n_sig = sum(1 for r in t4 if f(r["fdr"]) < 0.05)
print("  features with FDR<0.05: %d/13 (text states 12 of 13)" % n_sig)
checks.append(n_sig == 12)   # 2026-10-18 重跑后为 12/13；正文与此一致

print()
print("=" * 108)
print("B. Axis-feature correlations (results/immune_axis_correlation.tsv)")
print("=" * 108)
cor = load("immune_axis_correlation.tsv")
print("total pairs: %d" % len(cor))
allpos = all(f(r["rho"]) > 0 for r in cor)
print("all rho > 0:", allpos)
checks.append(allpos)
nut = [r for r in cor if r["axis"] == "nutrition_liver"]
print("nutrition_liver pairs: %d | zero-overlap pairs: %d"
      % (len(nut), sum(1 for r in nut if int(r["overlap_genes"]) == 0)))
checks.append(len(nut) == 13 and all(int(r["overlap_genes"]) == 0 for r in nut))
chk("nutrition rho min (NK cell activation)", min(f(r["rho"]) for r in nut), 0.069)
chk("nutrition rho max (mast cell activation)", max(f(r["rho"]) for r in nut), 0.362)
fdr_all = all(f(r["fdr"]) < 0.05 for r in nut)
print("all 13 nutrition pairs FDR<0.05:", fdr_all)
checks.append(fdr_all)
for target, claim in [("GOBP_MONOCYTE_ACTIVATION", 0.312), ("HALLMARK_INFLAMMATORY_RESPONSE", 0.209),
                      ("GOBP_MAST_CELL_ACTIVATION", 0.362), ("GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", 0.242)]:
    v = [f(r["rho"]) for r in nut if r["immune_feature"] == target][0]
    chk("nutrition rho vs %s" % target.replace("GOBP_", "").replace("HALLMARK_", "")[:26], v, claim)

print()
print("=" * 108)
print("C. Axis-transcript correlations (results/checkpoint_correlation.tsv)")
print("=" * 108)
ck = load("checkpoint_correlation.tsv")
print("total pairs: %d | transcripts: %d" % (len(ck), len({r["gene"] for r in ck})))
checks.append(len(ck) == 57)
cn = [r for r in ck if r["axis"] == "nutrition_liver"]
sig = [r for r in cn if f(r["p"]) < 0.05]
n_fdr = sum(1 for r in cn if f(r["fdr"]) < 0.05)
print("nutrition: n=%d | p<0.05: %d | FDR<0.05: %d (value quoted in the text) | rho range %.3f to %.3f"
      % (len(cn), len(sig), n_fdr, min(f(r["rho"]) for r in cn), max(f(r["rho"]) for r in cn)))
checks.append(len(cn) == 19 and len(sig) == 9 and n_fdr == 8)
chk("nutrition rho min across 19 transcripts", min(f(r["rho"]) for r in cn), -0.170)
chk("nutrition rho max across 19 transcripts", max(f(r["rho"]) for r in cn), 0.261)
for g, claim in [("CD274", -0.170), ("CTLA4", 0.063), ("LAG3", -0.032), ("PDCD1", 0.077),
                 ("TIGIT", 0.049), ("IFNG", -0.049)]:
    v = [f(r["rho"]) for r in cn if r["gene"] == g][0]
    chk("nutrition rho vs %s" % g, v, claim)

print()
print("=" * 108)
print("D. Univariate + multivariable Cox (LUAD)")
print("=" * 108)
# 基线：2026-10-18 R/08 重跑后重新校准（与 01_results.md 同步）
uni = {r["axis"]: r for r in load("axis_univariate_cox_luad.tsv")}
for a, hr, p in [("nutrition_liver", 1.090, 0.278), ("neutrophil", 0.963, 0.599), ("lymphocyte", 0.843, 0.025)]:
    chk("univariate HR %s" % a, f(uni[a]["HR"]), hr, tol=0.0015)
    chk("univariate P  %s" % a, f(uni[a]["p"]), p, tol=0.0015)

mv = {r["term"]: r for r in load("multivariate_cox_luad.tsv")}
for term, hr, lo, hi, p in [("risk_score", 15.568, 0.872, 277.839, 0.0619),
                            ("stage_latelate", 2.397, 1.749, 3.287, 5.6e-08),
                            ("age", 1.010, None, None, 0.204),
                            ("genderMALE", 1.015, None, None, 0.920)]:
    if term in mv:
        r = mv[term]
        chk("multivariable HR %s" % term, f(r["HR"]), hr)
        if lo is not None:
            chk("multivariable CI-lo %s" % term, f(r["lo"]), lo)
            chk("multivariable CI-hi %s" % term, f(r["hi"]), hi)
        if p < 1e-6:
            print("  %-52s %-28s %-24s %s" % ("multivariable P %s" % term, "%.1e" % f(r["p"]),
                                              "%.1e" % p, "OK" if abs(f(r["p"]) - p) / p < 0.05 else "**DIFF**"))

print()
print("=" * 108)
print("E. Model coefficients and lambda")
print("=" * 108)
mc = load("model_coef.tsv")
for r in mc:
    print("  %-16s coef=%+.4f  lambda_min=%.4g  s_used=%s"
          % (r["axis"], f(r["coef"]), f(r["lambda_min"]), r["s_used"]))
chk("lambda.min (full LUAD)", f(mc[0]["lambda_min"]), 0.0408, tol=0.00005, fmt="%.4f")
lym = [r for r in mc if r["axis"] == "lymphocyte"][0]
chk("retained coefficient (lymphocyte)", f(lym["coef"]), -0.053, tol=0.0005)
chk("zero-model fallback triggered", str(lym["s_used"] != "lambda.min"), "False")
chk("n non-zero coefficients", str(sum(1 for r in mc if abs(f(r["coef"])) > 1e-12)), "1")

print()
print("=" * 108)
print("F. Internal validation: single split + 100 splits")
print("=" * 108)
sp = load("internal_split_summary.tsv")[0]
print("single split:", dict(sp))
chk("single-split training C", f(sp["c_train"]), 0.587)
chk("single-split held-out C", f(sp["c_test"]), 0.494)
chk("single-split log-rank P", f(sp["p_test"]), 0.888)
sp_coefs = [f(x) for x in str(sp["coef"]).split("+")]
chk("single-split coef (nutrition_liver)", sp_coefs[0], 0.154)
chk("single-split coef (lymphocyte)", sp_coefs[1], -0.185)


def q7(x, p):
    """R quantile(type = 7): the method used for the manuscript percentiles."""
    n = len(x)
    h = (n - 1) * p
    lo = int(h)
    if lo + 1 >= n:
        return x[lo]
    return x[lo] + (h - lo) * (x[lo + 1] - x[lo])


cv = load("internal_repeated_cv.tsv")
ct = sorted(f(r["c_test"]) for r in cv)
gaps = [f(r["gap"]) for r in cv]
print("n splits: %d" % len(cv))
chk("median held-out C (100 splits)", st.median(ct), 0.542)
chk("mean held-out C", st.mean(ct), 0.541)
chk("2.5th percentile (type 7)", q7(ct, 0.025), 0.456, tol=0.002)
chk("97.5th percentile (type 7)", q7(ct, 0.975), 0.600, tol=0.002)
chk("median optimism", st.median(gaps), 0.017)
chk("mean optimism", st.mean(gaps), 0.015)
chk("n_sel == 1 in 92 of 100 splits", str(sum(1 for r in cv if int(r["n_sel"]) == 1)), "92")

aud = load("_audit_internal_selected.tsv")
from collections import Counter
cnt = Counter(r["selected"] for r in aud)
print("axis identity across splits:", dict(cnt))
nofall = sum(1 for r in aud if r["fellback"].strip().upper() == "TRUE")
print("splits with no axis passing screen (fellback): %d/100" % nofall)
lym_all = sum(v for k, v in cnt.items() if "lymphocyte" in k)
nut_all = sum(v for k, v in cnt.items() if "nutrition_liver" in k)
neu_all = sum(v for k, v in cnt.items() if "neutrophil" in k)
print("axis appears-in counts: lymphocyte %d, nutrition_liver %d, neutrophil %d"
      % (lym_all, nut_all, neu_all))
checks.append(cnt["lymphocyte"] == 82 and cnt["nutrition_liver"] == 8 and cnt["neutrophil"] == 2)
checks.append(lym_all == 90 and nut_all == 14 and neu_all == 4)
checks.append(nofall == 33)

print()
print("=" * 108)
print("G. Sensitivity analyses")
print("=" * 108)
gs = load("geneset_sensitivity.tsv")
print("schemes:", sorted({r["scheme"] for r in gs}), "| directions:", sorted({(r["train"], r["test"]) for r in gs}))
for r in gs:
    if r["train"] == "LUAD" and r["test"] == "LUSC":
        print("  LUAD->LUSC %-12s n_feat=%s n_sel=%s train=%.3f test=%.3f p=%s"
              % (r["scheme"], r["n_feature"], r["n_sel"], f(r["c_train"]), f(r["c_test"]), r["p_test"]))
    if r["train"] == "LUSC" and r["test"] == "LUAD":
        print("  LUSC->LUAD %-12s n_feat=%s n_sel=%s train=%.3f test=%.3f p=%s"
              % (r["scheme"], r["n_feature"], r["n_sel"], f(r["c_train"]), f(r["c_test"]), r["p_test"]))

vs = load("validation_summary.tsv")
for r in vs:
    print("  %-22s n=%s events=%s C=%.3f (%.3f-%.3f) logrank p=%.3f"
          % (r["cohort"], r["n"], r["events"], f(r["cindex"]), f(r["ci_lo"]), f(r["ci_hi"]), f(r["logrank_p"])))
lusc = [r for r in vs if "LUSC" in r["cohort"] and "train" not in r["cohort"]]
if lusc:
    chk("LUSC validation C-index", f(lusc[0]["cindex"]), 0.501)
    chk("LUSC validation log-rank P", f(lusc[0]["logrank_p"]), 0.751)

print()
print("=" * 108)
print("SUMMARY: %d/%d checks passed" % (sum(checks), len(checks)))
print("=" * 108)

# Analysis pipeline — run order, inputs and outputs

The analysis is written in **R** (R 4.6.1) and is driven by a single master script. A small number of
downstream tools are pure **Python** (3.13) and use only the standard library except where noted.

---

## 0. Master entry point

```bash
Rscript run_all_pathA.R                  # full pipeline
Rscript run_all_pathA.R 99               # synthetic self-test only
Rscript run_all_pathA.R 01 02 03 04      # selected steps only
```

`run_all_pathA.R` calls `setwd()` to the project root, then `source()`s `R/<step>_<name>.R` for each step.
**Full rebuild order: `01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11`.** Step `11` must run **before**
the article-table index is built, because the index reads `results/tables/TableS4{a,b,c,d}_*.csv`:

```bash
python tools/make_article_tables.py
```

> **Reproducing on another machine.** Every script contains one hard-coded project root. Run
> `python relocate.py` first — see README §5.

---

## 1. Configuration and shared library (not a step)

| File | Purpose |
|---|---|
| `R/00_setup_pathA.R` | Global configuration: all input/output paths (`PATHA`), analysis parameters (`PARS`: endpoint OS, training cohort TCGA-LUAD, validation cohort TCGA-LUSC, 10-fold LASSO, α = 1, 500 bootstrap replicates, seed 20260101), and the **three a-priori nutrition–inflammation axes** (`GENES`: 16 + 16 + 16 genes). Also creates the output directories. |
| `R/pathA_core.R` | Core function library: Toil probe-map reader, expression loader with `toil_fix()` (reverses the stored log2(FPKM + 0.001) encoding), Ensembl→symbol collapse (mean on the as-distributed scale), ssGSEA scoring via the GSVA 2.x `ssgseaParam` API, LASSO-Cox fitting, C-index evaluation, Kaplan–Meier helper. Written so that each function is pure and independently testable; deliberately does **not** depend on survminer / estimate / WGCNA / survAUC. |

---

## 2. Steps 01–13

| # | Script | Purpose | Key inputs | Key outputs |
|---|---|---|---|---|
| 01 | `R/01_load.R` | Read the Toil expression matrix + sample metadata + TCGA-CDR clinical table; build the LUAD (training) and LUSC (validation) sample-level survival tables from primary-tumour (`-01`) samples only | `data/tcga_RSEM_gene_fpkm.gz`, `work/expr_luad_lusc.tsv` (+ `.meta.tsv`), `assets/_gencode_probemap.gz`, `data/tcga_cdr_clinical.csv` | `work/pa_data.rds` (`expr_log` = log2(FPKM+1) symbol-level matrix; `surv_luad`; `surv_lusc`; `genesets`) |
| 02 | `R/02_score.R` | ssGSEA scoring of all primary-tumour samples against the three nutrition–inflammation axes | `work/pa_data.rds` | `work/pa_scores.rds` (sample × 3 axis scores) |
| 03 | `R/03_lasso_cox.R` | **LASSO-Cox modelling** in the LUAD training set (10-fold CV, α = 1); extract coefficients and the training median risk score | `work/pa_data.rds`, `work/pa_scores.rds` | `work/pa_model.rds`; `results/model_coef.tsv`; `results/axis_univariate_cox_luad.tsv` |
| 04 | `R/04_validate.R` | Three-level validation: LUAD internal (bootstrap C-index, KM, multivariable Cox) **and** LUSC cross-histology generalisation | `work/pa_*.rds` | `results/validation_summary.tsv`; `results/multivariate_cox_luad.tsv`; `figures/supplementary/FigureS1_km_luad.{png,pdf}` … `FigureS4_forest_luad.{png,pdf}` |
| 05 | `R/05_endpoint_sensitivity.R` | **Endpoint sensitivity**: OS / DSS / PFI run in parallel through the identical pipeline and gene sets, to test whether switching the endpoint rescues the signal | `work/pa_data.rds`, `work/pa_scores.rds` | `results/endpoint_sensitivity.tsv` |
| 06 | `R/06_geneset_sensitivity.R` | **Two variables tested against real data**: (A) feature scheme — S0 prior 3 axes (16 genes each) / S1 six real MSigDB pathways / S2 the 48 prior genes at single-gene level; (B) training direction — LUAD→LUSC **and** LUSC→LUAD. Endpoints OS and PFI each run once | `work/pa_*.rds`, `data/genesets/*.gmt` | `results/geneset_sensitivity.tsv` |
| 07 | `R/07_immune_microenvironment.R` | **The pivot of the paper**: nutrition–inflammation axes × tumour immune microenvironment. Pathway-level activity proxies from real GO:BP / HALLMARK entries (explicitly **not** a deconvolution estimate) | `work/pa_scores.rds`, expression matrix, local `.gmt` | `results/immune_axis_correlation.tsv`, `results/checkpoint_correlation.tsv`, `results/immune_group_diff.tsv`, `work/pa_scores_immune.rds`; *(intermediate heat map → `figures/_intermediate/`, not part of the paper)* |
| 08 | `R/08_internal_validation.R` | LUAD-internal train/test validation, adopted after Step 06 showed a cross-histology C-index ceiling of ≈0.55. Simple random 7:3 split (deliberately **not** event-stratified), feature selection and scaling done **inside the training set only**, then 100 repeated splits to report the distribution | `work/pa_*.rds` | `results/internal_split_summary.tsv`, `results/internal_repeated_cv.tsv`; *(intermediate KM plot → `figures/_intermediate/`)* |
| 09 | `R/09_tables.R` | Assemble the article-level tables from already-computed intermediate products — **no new analysis**. Uniform formatting rules (P values, one decimal for percentages, three decimals for rho) | `work/*.rds`, `results/*.tsv` | `results/tables/Table1_baseline.csv` … `Table5_*.csv`, `TableS1_endpoint_sensitivity.csv`, `TableS2_geneset_sensitivity.csv` |
| 10 | `R/10_figures_pathA.R` | Manuscript figures 2–5 (Fig. 1 is produced in Python, see §3). Layout compliant with the journal's artwork guidance: sans-serif 8–12 pt, no in-figure title or caption | `results/*.tsv`, `results/figure_data/*.tsv` | `figures/fig2_axis_immune_heatmap`, `fig3_nutrition_group_boxplot`, `fig4_nutrition_effect_forest`, `fig5_internal_validation` (`.png` + `.pdf`) |
| 11 | `R/11_tableS4.R` | Export Table S4 (parts a–d): the gene / pathway / transcript inventories moved out of the Methods text. **Reads only already-written source definitions and results — performs no new computation and changes no value**; every field traces back to `R/00_setup_pathA.R`, `R/06`, `R/07` or `results/*.tsv` | `R/00_setup_pathA.R`, `results/*.tsv` | `results/tables/TableS4a_axis_genes.csv`, `TableS4b_immune_features.csv`, `TableS4c_checkpoint_effector.csv`, `TableS4d_msigdb6_pathways.csv` |
| 12 | `R/12_random_panel_estimate.R` | Two robustness modules. (a) **Random-panel null models**: for each axis, 1,000 random 16-gene panels under an expression-matched null (twentieths of the expression-level distribution) and an unconstrained null, all 6,003 panels scored in the same ssGSEA call as the observed axes; compares mean \|Spearman rho\| across the 13 immune features, per-feature correlations, and single-axis OS C-indices. (b) **ESTIMATE purity adjustment**: matches precomputed ESTIMATE scores by 15-character barcode (1,010/1,011 samples) and recomputes every axis–feature association as a Spearman partial correlation controlling for `tumor.purity`. Requires `data/external/estimate/ESTIMATE_scores.txt` (DATA_SOURCES.md §6); skipped with a message otherwise. K=1,000; ~20 min | `work/expr_luad_lusc.tsv`, `work/pa_scores*.rds`, `data/external/estimate/*` | `results/random_panel_null_summary.tsv`, `random_panel_null_immune.tsv`, `random_panel_null_survival.tsv`, `random_panel_draws.tsv.gz`, `estimate_axis_correlation.tsv`, `estimate_partial_correlation.tsv`, `estimate_partial_summary.tsv`, `estimate_scores_matched.tsv` |
| 13 | `R/13_figures_random_estimate.R` | Main figures 6–7 from the step-12 outputs (no re-computation): Fig. 6 random-panel null distributions vs observed axes; Fig. 7 raw vs purity-adjusted correlations by gene-sharing class. A6-0 artwork standard (6.8 in wide, 600 dpi PNG + cairo PDF, no in-figure title) | `results/random_panel_draws.tsv.gz`, `results/random_panel_null_summary.tsv`, `results/estimate_partial_correlation.tsv` | `figures/fig6_random_null.{png,pdf}`, `figures/fig7_purity_adjust.{png,pdf}` |
| 99 | `R/99_selftest.R` | **Synthetic end-to-end self-test, independent of the TCGA data.** 150 samples × 60 genes, three gene sets (16/8/8), a positive survival-risk signal implanted on the neutrophil axis (HR ≈ 1.05 per unit score). Pass criterion: no error, C > 0.55, coefficient direction correct, consistent direction on the (synthetic) LUSC hold-out | — (synthetic) | console PASS/FAIL lines; exit status |

> **Note on step 06.** The header comment of `R/06_geneset_sensitivity.R` mentions a second output file,
> `results/axis_direction_double.tsv`. In the shipped version the training-direction (part B) results are
> written **into `results/geneset_sensitivity.tsv` together with part A**, and the standalone file is not
> produced. This is a stale comment. The script is shipped unmodified (see "Reproducibility" below), so the
> comment was left as-is; the authoritative list of produced files is this table.

---

## 3. Downstream tools (`tools/`)

### 3.1 Analysis-side

| Script | Purpose | Inputs | Outputs |
|---|---|---|---|
| `tools/subset_toil.py` | Stream-subsets the gzipped pan-cancer Toil matrix down to LUAD/LUSC (memory-friendly: line by line, keeping only the target sample columns) | `data/tcga_RSEM_gene_fpkm.gz`, `data/tcga_sample_barcodes.tsv` | `work/expr_luad_lusc.tsv` + `work/expr_luad_lusc.tsv.meta.tsv` |
| `tools/make_fig1_workflow.py` | **Figure 1** — study-design and analysis workflow diagram (matplotlib, 600 dpi PNG + embedded-font PDF) | numbers read from `results/tables/Table1_baseline.csv`, `results/internal_split_summary.tsv`, `results/internal_repeated_cv.tsv` | `figures/fig1_workflow.{png,pdf}` |
| `tools/make_article_tables.py` | Concatenates `results/tables/*.csv` into a single article-level table document | `results/tables/*.csv` | `results/tables/PathA_文章级表格.md` / `.html` |
| `tools/verify_msigdb_version.sh` | **Pins the MSigDB version** by byte-level md5 comparison against the Broad public release archive (read-only on the remote; writes only a cache under `data/_vercheck/`) | `data/genesets/*.gmt`, Broad archive | `data/_vercheck/msigdb_version_verdict.txt` |
| `tools/verify_table4.R` | Independently recomputes the nutrition-axis median-split comparison (article Table 4) to confirm reproducibility, and persists the **unrounded** values so every number in the text can be traced back; writes the same table under two names (`table4_…`/`table5_…`) for backward compatibility with earlier consumers | `work/pa_scores.rds`, `work/pa_scores_immune.rds` | `results/table4_nutrition_axis_raw.tsv` (= `results/table5_nutrition_axis_raw.tsv`) |
| `tools/rebuild_expr_subset.py` | Stream-rebuilds `work/expr_luad_lusc.tsv` from the raw Toil compendium **by barcode name** (row per gene, original string values, no float round-trip). Added in v1.1.0 as the verified repair path after a column-label misalignment was found and fixed in the v1.0.0-era subset | `data/tcga_RSEM_gene_fpkm.gz`, `data/tcga_sample_barcodes.tsv` | `work/expr_luad_lusc.tsv` (1,120 sample columns × 60,498 gene rows) |
| `tools/audit_methods_facts.R` | **Read-only fact audit** performed before the Methods section was written: every count and parameter quoted in Methods is measured back from the source products, never from memory | `work/pa_data.rds` | console report (also kept in `logs/audit_methods_facts.log`) |
| `tools/audit_internal_selected.R` | Audits which features were actually selected in the internal-validation repeats | `work/pa_*.rds` | `results/_audit_internal_selected.tsv` |
| `tools/audit_screen_pvalues.R` | Audits the screening P values quoted in Methods | `work/pa_data.rds` | console report |

### 3.2 Manuscript-side (auxiliary — require the `manuscript/` tree, which is **not** shipped here)

| Script | Purpose |
|---|---|
| `tools/epmc_gapfill.py` | Targeted literature gap-fill against Europe PMC (`count` / `search`); writes only records actually returned by the API — no field is invented or inferred |
| `tools/a5_refmeta.py`, `tools/a5_shortlist.py`, `tools/a5_renumber.py` | Reference metadata retrieval, shortlisting and renumbering for the Introduction/Discussion |
| `tools/a6_renumber.py` | Final manuscript-wide reference-renumbering pass with an idempotence guard and a word-count baseline assertion |
| `tools/a6_3_esm.py` | Rewrites supplementary-material citations to the journal's required **file-name** form (`ESM 1`–`ESM 9`), with a residual scan, a no-space-filename safety check and a `--check` dry run |
| `tools/md2html.py` | Shared Markdown→HTML renderer for the manuscript pages |

---

## 4. Notes on reproducibility

1. **Root path.** Every R script and several Python tools contain a hard-coded project root. Run
   `python relocate.py` once before re-running (README §5).
2. **Byte-level provenance.** All scripts are shipped **unmodified**, so the archived code is identical to
   the code that produced the published results. The single exception is the root literal, which
   `relocate.py` rewrites on the user's machine (and `--restore` reverts).
3. **Numbers are never hand-copied.** Every document that reports a number reads it directly from a result
   file. The visible text of the manuscript was checked against these files by an automated assertion
   suite (116 assertions).
4. **Zero-difference discipline.** Any re-run that changes only output paths or plotting style must leave
   the numerical content of every `results/*.tsv` byte-identical.
5. **Intermediate objects are not shipped.** `work/*.rds` (642 MB) are regenerated by steps 01–08 from the
   inputs listed in `DATA_SOURCES.md`.
6. **Raw data are not shipped** — see `DATA_SOURCES.md`.
7. **Integrity of the shipped files** is verifiable with `python verify_checksums.py` against
   `CHECKSUMS.sha256`.

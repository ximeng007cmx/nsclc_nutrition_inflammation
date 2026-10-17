# Analysis code — Nutrition–inflammation gene-expression axes associate with myeloid inflammatory activity but do not predict survival in non-small cell lung cancer

Version **v1.1.0** ｜ Date **2026-10-18**

Analysis code accompanying the manuscript:

> Nutrition–inflammation gene-expression axes associate with myeloid inflammatory activity but do not predict survival in non-small cell lung cancer: a TCGA transcriptomic analysis

This repository contains the **complete analysis pipeline** (data assembly → ssGSEA axis scoring → LASSO-Cox
modelling → internal and cross-histology validation → endpoint, feature-scheme and training-direction
sensitivity analyses → random-panel null calibration and ESTIMATE purity adjustment →
immune-microenvironment analysis → article-level tables and figures), together with
the **derived result tables** on which every reported number is based.

---

## 1. Study summary

- **Cohort**: TCGA primary lung tumours from the Toil pan-cancer RNA-seq compendium — **1,011 primary tumours**
  (TCGA-LUAD **n = 500** training; TCGA-LUSC **n = 491** cross-histology validation), 58,581 mapped gene symbols.
- **Exposure of interest**: three a-priori **nutrition–inflammation axis** signatures (16 genes each) built as
  proxies for the three components of the advanced lung cancer inflammation index (ALI) — a liver-derived
  nutritional axis, a neutrophil inflammatory axis, and a lymphocyte/cytotoxic-T axis. These are **prior
  definitions**, not a data-driven screen; the manuscript therefore references "nutrition–inflammation
  molecular signatures" and does **not** claim equivalence with the clinical ALI.
- **Primary question**: do these axes (a) track the tumour immune microenvironment, and (b) predict survival?
- **Headline findings** (all reproducible from `results/`):
  - The axes **do** associate with myeloid inflammatory activity (pathway-level proxies from real GO:BP and
    HALLMARK entries — explicitly **not** a deconvolution estimate).
  - They **do not** predict survival: internal hold-out C-index **0.494** (log-rank P = 0.888; single 7:3 split,
    n_test = 150) and a 100-split median test C-index of **0.542** (2.5th–97.5th percentile 0.456–0.600);
    cross-histology external C-index 0.501 (LUAD-fitted model applied to LUSC).
  - **Random-panel calibration**: the lymphocyte and neutrophil axes exceed 1,000 expression-matched random
    16-gene panels (empirical P = 0.001), whereas the nutrition axis lies within the null (P = 0.18).
  - **Purity adjustment**: partialling out the ESTIMATE tumour-purity estimate dissolves the lymphocyte-axis
    correlations with unrelated myeloid programmes (median rho 0.645 → −0.035) but leaves the nutrition-axis
    myeloid associations intact (median 0.207 → 0.057, 7 of 13 still FDR-significant).
  - The prognostic result is reported as an **honest negative**, not as a rescued model.
- **Design principles**: a-priori signatures; feature selection and scaling performed **inside** the training
  set only; sensitivity analyses for the endpoint (OS / DSS / PFI), the feature scheme (3 prior axes vs. 6 real
  MSigDB pathways vs. 48 single genes) and the training direction (LUAD→LUSC and LUSC→LUAD).

---

## 2. Repository layout

```
.
├── README.md                       this file
├── ANALYSIS_PIPELINE.md            run order, inputs and outputs of every script
├── DATA_SOURCES.md                 every external dataset, with accession, URL and fingerprint
├── CHECKSUMS.sha256                sha256 of every file shipped in this repository
├── verify_checksums.py             checks the tree against CHECKSUMS.sha256
├── requirements.txt                R and Python environment
├── relocate.py                     rewrite the hard-coded project root (see §5)
├── run_all_pathA.R                 master entry point (steps 01–13, 99)
├── LICENSE                         MIT
├── CITATION.cff                    citation metadata (Citation File Format 1.2.0)
├── .zenodo.json                    Zenodo deposition metadata for the archived release
├── R/                              the analysis steps
│   ├── 00_setup_pathA.R            configuration: paths, parameters, the three a-priori signatures
│   ├── pathA_core.R                shared function library (Toil handling, ssGSEA, LASSO-Cox, evaluation)
│   ├── 01_load.R … 13_figures_random_estimate.R   steps 01–13 (12–13 need the ESTIMATE scores, §4)
│   └── 99_selftest.R               synthetic end-to-end self-test (no TCGA data needed)
├── tools/                          downstream tools (Python + a few R audit scripts)
├── results/                        derived result tables cited in the manuscript
│   ├── tables/                     article-level Tables 1–5 and S1–S4
│   ├── figure_data/                the exact numbers behind Figure 4
│   └── refs_gapfill/               Europe PMC literature gap-fill records (manuscript auxiliary)
├── figures/                        Figures 1–7 (PNG + PDF)
│   └── supplementary/              Supplementary Figures S1–S4 (PNG + PDF)
├── logs/                           run logs + logs/sessionInfo.txt (full resolved R session)
└── assets/
    └── _gencode_probemap.gz        GENCODE gene-id → symbol annotation map (redistributed)
```

**Result files** (`results/`) are the derived numerical outputs behind the manuscript's tables and figures.
They are included so that the reported numbers can be checked **without** re-downloading the 720 MB
expression matrix or re-running the pipeline.

> **Not included**: the Toil expression matrix (720 MB), the TCGA-CDR clinical table, the sample barcode
> list, the MSigDB `.gmt` gene sets (all third-party — accession, URL and md5 in `DATA_SOURCES.md`) and the
> intermediate `work/*.rds` objects (642 MB, regenerated by steps 01–08).

---

## 3. Environment

- **R 4.6.1** with `data.table`, `GSVA` (2.6.6), `glmnet` (5.1), `survival` (3.8-6), `ggplot2`, `patchwork`, `pheatmap`
- **Python 3.13** with `matplotlib` (only `tools/make_fig1_workflow.py` needs it; the rest use the standard library)
- No command-line tool is required.
- The complete resolved environment is recorded in **`logs/sessionInfo.txt`**.

```bash
# R packages (Bioconductor for GSVA)
Rscript -e 'install.packages(c("data.table","glmnet","survival","ggplot2","patchwork","pheatmap"))'
Rscript -e 'if (!requireNamespace("BiocManager")) install.packages("BiocManager"); BiocManager::install("GSVA")'

# Python
python -m venv .venv
# Windows:  .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
```

---

## 4. How to run

```bash
# 0) point the scripts at your local checkout (see §5)
python relocate.py .

# 1) obtain the inputs listed in DATA_SOURCES.md, then place them as:
#      data/tcga_RSEM_gene_fpkm.gz
#      data/tcga_cdr_clinical.csv
#      data/tcga_sample_barcodes.tsv
#      data/genesets/h.all.symbols.gmt
#      data/genesets/c5.go.bp.symbols.gmt
#      data/external/estimate/ESTIMATE_scores.txt          (only needed for steps 12–13;
#      data/external/estimate/lung_*_RNAseqV2.txt           see DATA_SOURCES.md §6)
#      assets/_gencode_probemap.gz            (already shipped here)

# 2) subset the pan-cancer matrix to LUAD/LUSC
python tools/subset_toil.py data/tcga_RSEM_gene_fpkm.gz data/tcga_sample_barcodes.tsv \
    work/expr_luad_lusc.tsv --project TCGA-LUAD,TCGA-LUSC --suffix 01

# 3) run the pipeline (steps 01 -> 11; add 12 13 for the random-panel and purity analyses,
#    which run only if data/external/estimate/ESTIMATE_scores.txt is present)
Rscript run_all_pathA.R 01 02 03 04 05 06 07 08 09 10 11
Rscript run_all_pathA.R 12 13   # robustness modules (~20 min for step 12)

# 4) Figure 1 + the combined article-table document
python tools/make_fig1_workflow.py
python tools/make_article_tables.py

# (optional) check the archived files against CHECKSUMS.sha256
python verify_checksums.py
```

To verify the pipeline logic **without** any TCGA data (synthetic 150 × 60 data set with an implanted
signal), run `Rscript run_all_pathA.R 99`.

See **`ANALYSIS_PIPELINE.md`** for the full ordered list, inputs and outputs of every script.

---

## 5. ⚠ Relocating the project root

The scripts were written for a fixed working directory and contain one hard-coded root, e.g.

```r
setwd("D:/WorkBuddy_Y02/out_ali_directions/pathA")
```

This is preserved deliberately so that the archived code is **byte-identical** to the code that produced the
published results. To run the pipeline elsewhere, execute:

```bash
python relocate.py              # rewrites the root to this checkout's absolute path
python relocate.py --dry-run    # show what would change, write nothing
python relocate.py --restore    # revert to the original root
```

`relocate.py` rewrites only that single root literal and prints a per-file report. Re-running the analysis
with the relocated paths reproduces the archived `results/` files. Note that the hard-coded root appears in
**both** the R scripts and several Python tools, so run `relocate.py` before anything else.

---

## 6. Data integrity

```bash
python verify_checksums.py            # verify every shipped file against CHECKSUMS.sha256
python verify_checksums.py --list     # list the manifest without checking
```

`CHECKSUMS.sha256` covers every file committed to this repository. **`Bash tools/verify_msigdb_version.sh`**
additionally pins the MSigDB version against the Broad public archive (see `DATA_SOURCES.md` §5).

---

## 7. License

MIT — see `LICENSE`.

**Scope of the licence.** The MIT licence covers the **analysis code** in this repository only. The TCGA
expression matrix, the TCGA-CDR clinical table, the sample barcodes and the MSigDB gene sets are
third-party data released under their own terms; they are listed with their accessions and sources in
`DATA_SOURCES.md` and are **not** redistributed here. The GENCODE probe map in `assets/` **is** redistributed,
because the analysis's symbol mapping must stay reproducible even if the upstream hub copy is unavailable.

---

## 8. Citation

See `CITATION.cff`. If you use this code, please cite the accompanying manuscript and the archived Zenodo
release of this repository.

| Item | Value |
|---|---|
| GitHub repository | https://github.com/ximeng007cmx/nsclc_nutrition_inflammation |
| Concept DOI (always resolves to the newest archived version) | https://doi.org/10.5281/zenodo.23427752 |
| Version DOI **v1.1.0** | https://doi.org/10.5281/zenodo.23443637 |
| Version DOI **v1.0.0** | https://doi.org/10.5281/zenodo.23427753 |

Cite the **version DOI of the release you actually used**. The manuscript cites **v1.1.0**, the version that
reproduces every reported number (including the random-panel null calibration and the ESTIMATE purity
adjustment).

---

## 9. Contact

Mingxiang Chen — Department of Hepatobiliary, Vascular and Hernia Surgery, Shapingba Hospital Affiliated to
Chongqing University. ORCID [0000-0002-0997-0957](https://orcid.org/0000-0002-0997-0957).

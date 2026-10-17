# Data sources

All data analysed in this study are **public, de-identified, third-party resources**. No individual-level
data were newly generated. Raw input files are **not redistributed** in this repository; every accession,
download location and locally verified fingerprint is listed below so that the inputs can be re-obtained
and checked.

**All fingerprints in this file were measured on the local copies on 2026-10-16** (expression, barcodes,
clinical, probe map, gene sets; sizes with `stat`, hashes with `md5sum`) **and on 2026-10-17** (ESTIMATE
scores, §6). Files that are redistributed by this repository are additionally listed in
`CHECKSUMS.sha256`.

---

## 1. Gene-level expression — TCGA Pan-Cancer RNA-seq compendium (Toil hub)

| Item | Value |
|---|---|
| File (as distributed) | `tcga_RSEM_gene_fpkm.gz` |
| Distributing hub | UCSC Xena **Toil** hub |
| URL | `https://toil.xenahubs.net/download/tcga_RSEM_gene_fpkm.gz` |
| Size | **720,187,853 bytes** |
| Verified by | size + 1 MiB prefix hash (see note) |
| First 1 MiB md5 (local) | `7e91aec421fe8abd5b588d162f1e1dda` |
| First 1 MiB md5 (remote) | `7e91aec421fe8abd5b588d162f1e1dda` — **identical** |
| Remote `Last-Modified` | **Fri, 09 Apr 2021 20:04:52 GMT** |
| Encoding | log2(FPKM + 0.001), floor −9.9658 |
| Sample columns | 10,535 (matrix header has 10,536 fields incl. the gene column) |

**Strength of the provenance claim (stated honestly).** The local file and the remote object are
**byte-identical over the whole length and over the first 1 MiB**, which is what fixes the source URL. Because
the Xena hub does **not** publish a checksum manifest for this file (`md5sum.txt` returns HTTP 403
AccessDenied), a full-file md5 comparison against a published reference value is not possible; a full
720 MB transfer was therefore not performed. The claim is "**same source**", not "proved identical
byte-for-byte against a published checksum". If a stronger statement is needed, download the full file and
compare local `md5sum` output against the remote hash obtained at that time.

**Reference.** Vivian J, Rao AA, Nothaft FA, et al. *Nat Biotechnol* 2017;35(4):314–316.
**PMID 28398314**; DOI [10.1038/nbt.3772](https://doi.org/10.1038/nbt.3772).
Goldman MJ, Craft B, Hastie M, et al. *Nat Biotechnol* 2020;38(6):675–678.
**PMID 32444850**; DOI [10.1038/s41587-020-0546-8](https://doi.org/10.1038/s41587-020-0546-8).

---

## 2. Sample barcodes

| Item | Value |
|---|---|
| Local file | `data/tcga_sample_barcodes.tsv` (not redistributed) |
| Size | 98,359 bytes |
| md5 | `365e13d2fad7ca59477dbebd43ed378f` |
| Source | **Genomic Data Commons (GDC) API** |
| Processing | reduced to the first four barcode fields, fourth field truncated to its two-character sample-type code (`TCGA-44-3918-01A` → `TCGA-44-3918-01`); duplicate keys collapsed to the first occurrence |

---

## 3. Clinical annotation — TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR)

| Item | Value |
|---|---|
| Local file | `data/tcga_cdr_clinical.csv` (not redistributed) |
| Size | 249,118 bytes |
| md5 | `2169c15c8cf866a82c3887d06298227b` |
| Provides | curated OS, DSS and PFI times and event indicators |

**Reference.** Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data
resource to drive high-quality survival outcome analytics. *Cell* 2018;173(2):400–416.e11.
**PMID 29625055**; DOI [10.1016/j.cell.2018.02.052](https://doi.org/10.1016/j.cell.2018.02.052).

---

## 4. GENCODE probe map (Ensembl gene id → HGNC symbol)

| Item | Value |
|---|---|
| Local file | `assets/_gencode_probemap.gz` — **redistributed here** (see `CHECKSUMS.sha256`) |
| Size | 3,244,244 bytes |
| md5 | `c9f7e859d17d4c9102c07b98ffbd753d` |
| Rows | 60,498 data rows (60,499 lines incl. header) |
| Columns | `id  gene  chrom  chromStart  chromEnd  strand` (6) |
| Distinct symbols | 58,581 (the value reported in Methods §2.2) |

The file is the six-column gene annotation map in the format used by the Toil hub
(`gencode.v23.annotation.gene.probemap`). **Provenance is reported as not byte-verified.** A byte-level
comparison against the hub copy was attempted on 2026-10-16 but the remote object is no longer publicly
readable: `https://toil.xenahubs.net/download/gencode.v23.annotation.gene.probemap` redirects to S3 and
returns **403 AccessDenied**, and the hub's `md5sum.txt` manifest also returns 403. The file is therefore
identified only by its own fingerprint above, and the symbol count it yields (58,581) is the value used
throughout the manuscript.

> Because this small annotation file is redistributed, the symbol mapping used by the analysis is fully
> reproducible even if the hub copy stays unavailable. Nothing in the analysis depends on gene coordinates
> or on the specific GENCODE build beyond the id→symbol mapping.

---

## 5. Gene sets — MSigDB

| Local file | Collection | Size | md5 | Version |
|---|---|---|---|---|
| `data/genesets/h.all.symbols.gmt` | HALLMARK (H) | 48,690 B | `d3606c325b9862bf3ea7991c2c2bf7d7` | **v2024.1.Hs** |
| `data/genesets/c5.go.bp.symbols.gmt` | C5 GO:BP | 4,940,222 B | `fd535c44825f9be3ea7c5d31e3d45f1d` | **v2024.1.Hs** |

**The version is confirmed, not inferred.** Both files were compared **byte-for-byte (md5)** against the
Broad Institute's *public* release archive, which requires no login:

```
https://data.broadinstitute.org/gsea-msigdb/msigdb/release/<version>/<collection>.v<version>.symbols.gmt
```

Six candidate versions were tested (7.5.1, 2023.1.Hs, 2023.2.Hs, 2024.1.Hs, 2025.1.Hs, 2026.1.Hs); each
local file matched **exactly one** of them, namely **2024.1.Hs**. Note that the 2023.2.Hs and 2024.1.Hs
HALLMARK files have the *same size* (48,690 B) but different md5 — size alone is not sufficient to pin a
version. The reproducible command is `bash tools/verify_msigdb_version.sh`.

**References.** Subramanian A, Tamayo P, Mootha VK, et al. *PNAS* 2005;102(43):15545–15550.
**PMID 16199517**; DOI [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102).
Liberzon A, Birger C, Thorvaldsdóttir H, et al. *Cell Syst* 2015;1(6):417–425.
**PMID 26771021**; DOI [10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004).

---

## 6. Tumour purity — precomputed ESTIMATE scores

| Local file | Content | Size | md5 |
|---|---|---|---|
| `data/external/estimate/ESTIMATE_scores.txt` | master table (all TCGA diseases): ID, Platform, Disease, ImmuneScore, StromalScore, ESTIMATEScore, tumor.purity; 13,117 data rows | 1,071,120 B | `f8607b208df23eaf491fa1a59b7883a2` |
| `data/external/estimate/lung_adenocarcinoma_RNAseqV2.txt` | TCGA-LUAD RNAseqV2 subset (517 data rows) | 19,985 B | `0c68a0e625b798a1a1e1055bdc26a245` |
| `data/external/estimate/lung_squamous_cell_carcinoma_RNAseqV2.txt` | TCGA-LUSC RNAseqV2 subset (501 data rows) | 19,605 B | `896ce7ca34f3b29f75b041992a7486a6` |

| Item | Value |
|---|---|
| Source | ESTIMATE authors' distribution site (Yoshihara et al. 2013) |
| URL pattern | `https://bioinformatics.mdanderson.org/estimate/tables/<disease>_<platform>.txt` and `https://bioinformatics.mdanderson.org/estimate/rawdata/ESTIMATE_scores.txt` |
| Matching | 15-character barcode prefix; 1,010 of the 1,011 primary tumours matched |
| Used for | Spearman partial correlations controlling for `tumor.purity` (`R/12_random_panel_estimate.R`, §4) |

These are third-party precomputed scores and are **not redistributed**. The per-disease tables were used
only to confirm row-level consistency with the master table (max absolute score difference 0 across all
matched rows); all analysis numbers come from the master table.

**Reference.** Yoshihara K, Shahmoradgoli M, Martínez E, et al. Inferring tumour purity and stromal and
immune cell admixture from expression data. *Nat Commun* 2013;4:2612.
**PMID 24113773**; DOI [10.1038/ncomms3612](https://doi.org/10.1038/ncomms3612).

---

## 7. Analysis software

| Tool | Version used | Used for |
|---|---|---|
| R | **4.6.1** (2026-06-24) | entire pipeline |
| GSVA | **2.6.6** | ssGSEA single-sample scoring (`ssgseaParam` + `gsva`) |
| glmnet | **5.1** | LASSO-Cox model (ties = "efron") |
| survival | **3.8-6** | Cox models, Kaplan–Meier, concordance (Harrell C) |
| data.table | **1.18.6.1** | I/O and reshaping |
| ggplot2 | **4.0.3** | figures |
| patchwork | **1.3.2** | figure composition |
| pheatmap | **1.0.13** | heatmaps |
| Python | **3.13.12** | downstream tools (`tools/*.py`) |
| matplotlib | ≥ 3.8 | Figure 1 flowchart only |

The complete resolved environment is recorded verbatim in `logs/sessionInfo.txt`.

**Method references.** Barbie DA, et al. *Nature* 2009;462(7269):108–112. **PMID 19847166**;
DOI [10.1038/nature08460](https://doi.org/10.1038/nature08460) (ssGSEA).
Hänzelmann S, Castelo R, Guinney J. *BMC Bioinformatics* 2013;14:7. **PMID 23323831**;
DOI [10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7) (GSVA).
Friedman J, Hastie T, Tibshirani R. *J Stat Softw* 2010;33(1):1–22. **PMID 20808728**;
DOI [10.18637/jss.v033.i01](https://doi.org/10.18637/jss.v033.i01) (glmnet).
Harrell FE, Lee KL, Mark DB. *Stat Med* 1996;15(4):361–387. **PMID 8668867**;
DOI [10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4](https://doi.org/10.1002/(SICI)1097-0258(19960229)15:4%3C361::AID-SIM168%3E3.0.CO;2-4)
(concordance index).

---

## 8. Data not redistributed here

| Item | Size | Reason |
|---|---|---|
| `tcga_RSEM_gene_fpkm.gz` | 720,187,853 B | third-party data; publicly downloadable from the URL above |
| `data/tcga_cdr_clinical.csv` | 249,118 B | third-party data (TCGA-CDR); download from the source publication's supplement |
| `data/tcga_sample_barcodes.tsv` | 98,359 B | derived from the GDC API (self-service query) |
| `data/external/estimate/*.txt` (three files) | 1.1 MB total | third-party precomputed ESTIMATE scores; download from the URLs in §6 |
| Intermediate `.rds` objects (`work/`, 642 MB) | 642 MB | regenerated by the pipeline from the inputs above |

The MSigDB `.gmt` gene sets are **also third-party**; they are not redistributed either, because the version
is pinned by md5 and the files can be re-downloaded deterministically from the Broad public archive
(§5). If a reviewer requires the exact gene-set files, they can be reconstructed byte-for-byte from the
URL pattern in §5.

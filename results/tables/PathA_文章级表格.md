# Path A — 文章级表格（营养-炎症分子特征 × 非小细胞肺癌免疫微环境）

> 生成日期：2026-10-16  
> 数据来源：UCSC Xena Toil（TCGA RSEM FPKM，经 log2(FPKM+1) 转换）+ TCGA-CDR 临床表 + MSigDB v2024.1.Hs  
> 分析样本：TCGA-LUAD 500 例（182 事件）、TCGA-LUSC 491 例（211 事件），均为原发肿瘤  
> 表述约束：本方案为 ALI 三组分的**分子代理**，不得表述为「等价于 ALI」；
> 免疫特征为**通路级活性代理**，不得表述为 CIBERSORT/ESTIMATE 去卷积结果。

## Table 1. Baseline characteristics of the study cohorts

TCGA-LUAD (training cohort) and TCGA-LUSC (cross-histology validation cohort). Clinical data from TCGA-CDR (Liu J et al., Cell 2018; PMID 29625055). Signature scores are ssGSEA values of the three nutrition–inflammation axes.

| Characteristic | LUAD | LUSC |
|---|---|---|
| Age, median (IQR) | 66 (59–72) | 68 (62–73) |
| Deaths, n (%) | 182 (36.4) | 211 (43.0) |
| Follow-up of censored, median (IQR), days | 670 (469–1129) | 749 (395–1367) |
| Male, n (%) | 230 (46.0) | 364 (74.1) |
| OS time, median (IQR), days | 656 (423–1127) | 660 (323–1264) |
| Patients, n | 500 | 491 |
| Stage early (I–II), n (%) | 388 (77.6) | 397 (80.9) |
| Stage late (III–IV), n (%) | 104 (20.8) | 90 (18.3) |
| Stage other/unclassified, n (%) | 8 (1.6) | 4 (0.8) |
| Stage unknown, n (%) | 0 (0.0) | 0 (0.0) |
| lymphocyte score, median (IQR) | 0.561 (0.482–0.633) | 0.558 (0.474–0.642) |
| neutrophil score, median (IQR) | 0.329 (0.257–0.403) | 0.326 (0.256–0.391) |
| nutrition_liver score, median (IQR) | 0.183 (0.108–0.302) | 0.202 (0.116–0.307) |


## Table 2. Correlation between nutrition–inflammation axes and immune features

Spearman correlations across 1,011 primary tumour samples. Immune features are MSigDB pathway-level activity proxies (GO Biological Process and Hallmark; see Table S4b for set sizes and checksums), NOT deconvolution estimates. 'No overlapping genes' indicates the axis and the immune feature share no genes, so the correlation is not definitional. FDR = Benjamini–Hochberg.

| Nutrition–inflammation axis | Immune feature (MSigDB) | Set size | Overlapping genes | Spearman rho | P | FDR | Type |
|---|---|---|---|---|---|---|---|
| lymphocyte | HALLMARK_ALLOGRAFT_REJECTION | 195 | 11 | 0.922 | <1e-300 | <1e-300 | Overlapping genes |
| lymphocyte | GOBP_T_CELL_ACTIVATION | 566 | 10 | 0.905 | <1e-300 | <1e-300 | Overlapping genes |
| lymphocyte | GOBP_REGULATORY_T_CELL_DIFFERENTIATION | 40 | 1 | 0.896 | <1e-300 | <1e-300 | Overlapping genes |
| lymphocyte | GOBP_NATURAL_KILLER_CELL_ACTIVATION | 97 | 0 | 0.865 | <1e-300 | <1e-300 | No overlapping genes |
| lymphocyte | GOBP_B_CELL_ACTIVATION | 269 | 2 | 0.847 | <1e-300 | <1e-300 | Overlapping genes |
| lymphocyte | GOBP_MYELOID_DENDRITIC_CELL_ACTIVATION | 29 | 0 | 0.840 | <1e-300 | <1e-300 | No overlapping genes |
| lymphocyte | HALLMARK_INTERFERON_GAMMA_RESPONSE | 196 | 2 | 0.831 | <1e-300 | <1e-300 | Overlapping genes |
| lymphocyte | GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION | 113 | 2 | 0.812 | <1e-300 | <1e-300 | Overlapping genes |
| lymphocyte | HALLMARK_INFLAMMATORY_RESPONSE | 199 | 5 | 0.704 | <1e-300 | <1e-300 | Overlapping genes |
| lymphocyte | GOBP_MAST_CELL_ACTIVATION | 70 | 0 | 0.649 | <1e-300 | <1e-300 | No overlapping genes |
| lymphocyte | GOBP_MACROPHAGE_ACTIVATION | 109 | 0 | 0.646 | <1e-300 | <1e-300 | No overlapping genes |
| lymphocyte | GOBP_NEUTROPHIL_DEGRANULATION | 14 | 0 | 0.590 | <1e-300 | <1e-300 | No overlapping genes |
| lymphocyte | GOBP_MONOCYTE_ACTIVATION | 12 | 0 | 0.421 | <1e-300 | <1e-300 | No overlapping genes |
| neutrophil | HALLMARK_INFLAMMATORY_RESPONSE | 199 | 1 | 0.629 | <1e-300 | <1e-300 | Overlapping genes |
| neutrophil | GOBP_MAST_CELL_ACTIVATION | 70 | 1 | 0.624 | <1e-300 | <1e-300 | Overlapping genes |
| neutrophil | GOBP_NEUTROPHIL_DEGRANULATION | 14 | 0 | 0.611 | <1e-300 | <1e-300 | No overlapping genes |
| neutrophil | GOBP_T_CELL_ACTIVATION | 566 | 1 | 0.598 | <1e-300 | <1e-300 | Overlapping genes |
| neutrophil | GOBP_B_CELL_ACTIVATION | 269 | 0 | 0.572 | <1e-300 | <1e-300 | No overlapping genes |
| neutrophil | GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION | 113 | 0 | 0.566 | <1e-300 | <1e-300 | No overlapping genes |
| neutrophil | HALLMARK_ALLOGRAFT_REJECTION | 195 | 2 | 0.560 | <1e-300 | <1e-300 | Overlapping genes |
| neutrophil | GOBP_MYELOID_DENDRITIC_CELL_ACTIVATION | 29 | 0 | 0.524 | <1e-300 | <1e-300 | No overlapping genes |
| neutrophil | GOBP_NATURAL_KILLER_CELL_ACTIVATION | 97 | 0 | 0.520 | <1e-300 | <1e-300 | No overlapping genes |
| neutrophil | GOBP_MACROPHAGE_ACTIVATION | 109 | 0 | 0.516 | <1e-300 | <1e-300 | No overlapping genes |
| neutrophil | HALLMARK_INTERFERON_GAMMA_RESPONSE | 196 | 1 | 0.497 | <1e-300 | <1e-300 | Overlapping genes |
| neutrophil | GOBP_REGULATORY_T_CELL_DIFFERENTIATION | 40 | 0 | 0.456 | <1e-300 | <1e-300 | No overlapping genes |
| neutrophil | GOBP_MONOCYTE_ACTIVATION | 12 | 0 | 0.335 | <1e-300 | <1e-300 | No overlapping genes |
| nutrition_liver | GOBP_MACROPHAGE_ACTIVATION | 109 | 0 | 0.276 | 5.3e-19 | 7.6e-19 | No overlapping genes |
| nutrition_liver | GOBP_MONOCYTE_ACTIVATION | 12 | 0 | 0.255 | 2.6e-16 | 3.6e-16 | No overlapping genes |
| nutrition_liver | HALLMARK_INFLAMMATORY_RESPONSE | 199 | 0 | 0.250 | 8.1e-16 | 1.1e-15 | No overlapping genes |
| nutrition_liver | GOBP_MAST_CELL_ACTIVATION | 70 | 0 | 0.226 | 4.6e-13 | 5.9e-13 | No overlapping genes |
| nutrition_liver | GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION | 113 | 0 | 0.197 | 3.2e-10 | 4.0e-10 | No overlapping genes |
| nutrition_liver | HALLMARK_ALLOGRAFT_REJECTION | 195 | 0 | 0.176 | 1.7e-08 | 2.1e-08 | No overlapping genes |
| nutrition_liver | GOBP_T_CELL_ACTIVATION | 566 | 0 | 0.165 | 1.4e-07 | 1.7e-07 | No overlapping genes |
| nutrition_liver | GOBP_MYELOID_DENDRITIC_CELL_ACTIVATION | 29 | 0 | 0.152 | 1.2e-06 | 1.4e-06 | No overlapping genes |
| nutrition_liver | HALLMARK_INTERFERON_GAMMA_RESPONSE | 196 | 0 | 0.144 | 4.6e-06 | 5.2e-06 | No overlapping genes |
| nutrition_liver | GOBP_B_CELL_ACTIVATION | 269 | 0 | 0.114 | 2.7e-04 | 3.0e-04 | No overlapping genes |
| nutrition_liver | GOBP_NATURAL_KILLER_CELL_ACTIVATION | 97 | 0 | 0.108 | 5.9e-04 | 6.2e-04 | No overlapping genes |
| nutrition_liver | GOBP_REGULATORY_T_CELL_DIFFERENTIATION | 40 | 0 | 0.100 | 0.002 | 0.002 | No overlapping genes |
| nutrition_liver | GOBP_NEUTROPHIL_DEGRANULATION | 14 | 0 | 0.099 | 0.002 | 0.002 | No overlapping genes |


## Table 3. Correlation between nutrition–inflammation axes and immune checkpoint / effector genes

Spearman correlations with single-gene expression (log2(FPKM+1)). Checkpoint panel: PDCD1, CD274, PDCD1LG2, CTLA4, LAG3, HAVCR2, TIGIT, IDO1, SIGLEC15. Effector panel: IFNG, GZMA, GZMB, PRF1, CXCL9, CXCL10, CXCL11, HLA-DRA, STAT1, CCR5.

| gene | panel | nutrition_liver (rho) | neutrophil (rho) | lymphocyte (rho) | nutrition_liver P | neutrophil P | lymphocyte P |
|---|---|---|---|---|---|---|---|
| CCR5 | Effector / TIS-like | 0.142 | 0.335 | 0.807 | 6.3e-06 | 6.1e-28 | 1.8e-233 |
| CXCL10 | Effector / TIS-like | 0.059 | 0.252 | 0.733 | 0.062 | 4.1e-16 | 3.3e-171 |
| CXCL11 | Effector / TIS-like | 0.031 | 0.245 | 0.675 | 0.319 | 3.1e-15 | 2.4e-135 |
| CXCL9 | Effector / TIS-like | 0.087 | 0.261 | 0.835 | 0.005 | 3.8e-17 | 1.6e-264 |
| GZMA | Effector / TIS-like | 0.156 | 0.379 | 0.914 | 6.7e-07 | 7.6e-36 | <1e-300 |
| GZMB | Effector / TIS-like | 0.070 | 0.454 | 0.838 | 0.026 | 1.5e-52 | 3.5e-267 |
| HLA-DRA | Effector / TIS-like | 0.156 | 0.370 | 0.731 | 6.2e-07 | 4.2e-34 | 6.1e-170 |
| IFNG | Effector / TIS-like | 0.033 | 0.303 | 0.814 | 0.297 | 6.6e-23 | 1.4e-240 |
| PRF1 | Effector / TIS-like | 0.116 | 0.404 | 0.870 | 2.2e-04 | 7.0e-41 | <1e-300 |
| STAT1 | Effector / TIS-like | -0.014 | 0.213 | 0.587 | 0.660 | 7.7e-12 | 1.1e-94 |
| CD274 | Immune checkpoint | -0.026 | 0.284 | 0.565 | 0.409 | 3.6e-20 | 2.3e-86 |
| CTLA4 | Immune checkpoint | 0.038 | 0.423 | 0.800 | 0.233 | 4.0e-45 | 3.6e-226 |
| HAVCR2 | Immune checkpoint | 0.142 | 0.255 | 0.522 | 5.4e-06 | 1.6e-16 | 9.3e-72 |
| IDO1 | Immune checkpoint | 0.098 | 0.297 | 0.732 | 0.002 | 4.5e-22 | 4.0e-170 |
| LAG3 | Immune checkpoint | 0.027 | 0.243 | 0.766 | 0.397 | 4.5e-15 | 2.0e-195 |
| PDCD1 | Immune checkpoint | 0.052 | 0.345 | 0.866 | 0.096 | 1.2e-29 | <1e-300 |
| PDCD1LG2 | Immune checkpoint | 0.066 | 0.323 | 0.666 | 0.037 | 6.1e-26 | 1.6e-130 |
| SIGLEC15 | Immune checkpoint | 0.059 | 0.157 | 0.185 | 0.059 | 5.6e-07 | 3.0e-09 |
| TIGIT | Immune checkpoint | 0.057 | 0.406 | 0.899 | 0.071 | 2.6e-41 | <1e-300 |


## Table 4. Internal validation of the prognostic score (negative result)

LUAD split 7:3 at RANDOM (sample.int(); NOT stratified by event status — observed event rates 38.3% training vs 32.0% held-out, against 36.4% expected under stratification); feature selection (univariate Cox P<0.10) and standardisation performed within the training set only. The held-out C-index is at chance level, indicating the three a-priori axes carry no prognostic discrimination in LUAD; the previously reported C-index of 0.541 was an over-optimistic value from fitting on the full dataset.

| Analysis | Value | Note |
|---|---|---|
| Single 7:3 split — training set | 0.536 | Features selected within training set only |
| Single 7:3 split — held-out test set | 0.488 (log-rank P = 0.822) | Selected by univariate Cox (P<0.10) then LASSO-Cox |
| Repeated 100 splits — test set, median [2.5th–97.5th] | 0.501 [0.438–0.560] | Same procedure re-run on 100 random splits |
| Mean optimism (training − test C) | 0.031 | Shows the degree of over-optimism in internally fitted C-index |


## Table 5. Immune features by nutrition-axis group (independent of definitional overlap)

Samples dichotomised at the median of the nutrition_liver axis, which shares NO genes with any immune feature. Higher nutrition-axis expression is associated with higher activity of myeloid/monocyte/macrophage and inflammatory pathways. Wilcoxon rank-sum test; FDR = Benjamini–Hochberg.

| immune_feature | Median (nutrition-high) | Median (nutrition-low) | Difference | P | FDR |
|---|---|---|---|---|---|
| GOBP_MONOCYTE_ACTIVATION | 1.160 | 1.106 | 0.054 | 6.7e-16 | 2.2e-15 |
| HALLMARK_INFLAMMATORY_RESPONSE | 1.081 | 1.032 | 0.050 | 4.4e-16 | 1.9e-15 |
| HALLMARK_ALLOGRAFT_REJECTION | 1.123 | 1.083 | 0.040 | 2.0e-09 | 4.2e-09 |
| GOBP_MAST_CELL_ACTIVATION | 0.945 | 0.905 | 0.040 | 4.4e-16 | 1.9e-15 |
| GOBP_MACROPHAGE_ACTIVATION | 1.003 | 0.968 | 0.035 | <1e-300 | <1e-300 |
| GOBP_MYELOID_DENDRITIC_CELL_ACTIVATION | 1.053 | 1.018 | 0.034 | 5.0e-08 | 8.1e-08 |
| GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION | 1.250 | 1.223 | 0.026 | 1.4e-10 | 3.7e-10 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 1.305 | 1.280 | 0.026 | 6.6e-07 | 9.5e-07 |
| GOBP_REGULATORY_T_CELL_DIFFERENTIATION | 0.978 | 0.952 | 0.025 | 1.1e-04 | 1.1e-04 |
| GOBP_NATURAL_KILLER_CELL_ACTIVATION | 0.753 | 0.728 | 0.025 | 1.9e-05 | 2.2e-05 |
| GOBP_T_CELL_ACTIVATION | 1.072 | 1.049 | 0.023 | 3.9e-09 | 7.2e-09 |
| GOBP_B_CELL_ACTIVATION | 1.009 | 0.988 | 0.021 | 2.9e-06 | 3.8e-06 |
| GOBP_NEUTROPHIL_DEGRANULATION | 1.215 | 1.206 | 0.009 | 8.0e-05 | 8.7e-05 |


## Table S1. Sensitivity analysis across endpoints (OS / DSS / PFI)

Same pipeline re-run with each endpoint. LUAD is essentially flat across endpoints (C-index 0.535–0.542); LUSC shows a significant association for DSS. Endpoints therefore do not rescue the LUAD training signal.

| Endpoint | Cohort | n | Events | C-index | P | Axes selected | Note |
|---|---|---|---|---|---|---|---|
| OS | LUAD | 500 | 182 | 0.541 | 0.100 | 3 | strongest single axis P=0.353; s_used=0.000504 |
| OS | LUSC | 491 | 211 | 0.544 | 0.070 | 1 | strongest single axis P=0.0462; s_used=lambda.min |
| DSS | LUAD | 465 | 112 | 0.535 | 0.871 | 3 | strongest single axis P=0.223; s_used=0.0003842 |
| DSS | LUSC | 439 | 86 | 0.569 | 0.003 | 1 | strongest single axis P=0.0335; s_used=lambda.min |
| PFI | LUAD | 500 | 206 | 0.542 | 0.055 | 3 | strongest single axis P=0.0379; s_used=lambda.min |
| PFI | LUSC | 491 | 143 | 0.561 | 0.080 | 3 | strongest single axis P=0.187; s_used=0.0004298 |


## Table S2. Sensitivity analysis across feature schemes and training directions

S0 = three a-priori axes (16 genes each); S1 = six MSigDB pathways (1,325 genes); S2 = 48 prior genes entered directly into LASSO-Cox. Internal C-index rises with feature number while external C-index falls — a direct demonstration of overfitting. Cross-histology validation ceiling ≈ 0.55.

| Endpoint | Training | Validation | Scheme | n features | n selected | C train | C test | P |
|---|---|---|---|---|---|---|---|---|
| OS | LUAD | LUSC | S0_prior16 | 3 | 3 | 0.541 | 0.550 | 0.560 |
| OS | LUAD | LUSC | S1_msigdb6 | 6 | 5 | 0.583 | 0.518 | 0.677 |
| OS | LUAD | LUSC | S2_genes48 | 48 | 47 | 0.659 | 0.516 | 0.728 |
| OS | LUSC | LUAD | S0_prior16 | 3 | 1 | 0.544 | 0.525 | 0.797 |
| OS | LUSC | LUAD | S1_msigdb6 | 6 | 6 | 0.531 | 0.514 | 0.799 |
| OS | LUSC | LUAD | S2_genes48 | 48 | 1 | 0.557 | 0.477 | 0.869 |
| PFI | LUAD | LUSC | S0_prior16 | 3 | 3 | 0.542 | 0.550 | 0.130 |
| PFI | LUAD | LUSC | S1_msigdb6 | 6 | 3 | 0.541 | 0.534 | 0.268 |
| PFI | LUAD | LUSC | S2_genes48 | 48 | 1 | 0.545 | 0.530 | 0.474 |
| PFI | LUSC | LUAD | S0_prior16 | 3 | 3 | 0.561 | 0.522 | 0.781 |
| PFI | LUSC | LUAD | S1_msigdb6 | 6 | 6 | 0.516 | 0.533 | 0.204 |
| PFI | LUSC | LUAD | S2_genes48 | 48 | 46 | 0.652 | 0.541 | 0.512 |


## Table S3. Immune features by composite risk score (definitional consistency check)

CAUTION: the composite score contains the lymphocyte axis with a negative coefficient, so lower immune activity in the high-score group is partly definitional. This table is provided for completeness only; Table 5 is the non-circular analysis.

| immune_feature | Median (high score) | Median (low score) | Difference | P | FDR |
|---|---|---|---|---|---|
| GOBP_NEUTROPHIL_DEGRANULATION | 1.196 | 1.218 | -0.022 | 3.0e-06 | 3.0e-06 |
| GOBP_MAST_CELL_ACTIVATION | 0.906 | 0.942 | -0.035 | 1.4e-12 | 1.5e-12 |
| GOBP_MACROPHAGE_ACTIVATION | 0.963 | 1.009 | -0.046 | 8.3e-24 | 1.2e-23 |
| GOBP_B_CELL_ACTIVATION | 0.974 | 1.020 | -0.046 | 1.7e-34 | 3.1e-34 |
| GOBP_MONOCYTE_ACTIVATION | 1.106 | 1.160 | -0.054 | 2.4e-16 | 2.9e-16 |
| GOBP_T_CELL_ACTIVATION | 1.031 | 1.086 | -0.056 | 3.6e-41 | 1.6e-40 |
| HALLMARK_INFLAMMATORY_RESPONSE | 1.026 | 1.086 | -0.060 | 3.2e-19 | 4.1e-19 |
| GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION | 1.204 | 1.264 | -0.060 | 5.1e-34 | 8.3e-34 |
| GOBP_NATURAL_KILLER_CELL_ACTIVATION | 0.711 | 0.780 | -0.069 | 2.7e-38 | 5.9e-38 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 1.264 | 1.334 | -0.069 | 8.8e-39 | 2.3e-38 |
| GOBP_MYELOID_DENDRITIC_CELL_ACTIVATION | 0.991 | 1.078 | -0.087 | 2.7e-39 | 8.6e-39 |
| GOBP_REGULATORY_T_CELL_DIFFERENTIATION | 0.923 | 1.011 | -0.088 | 6.6e-52 | 8.6e-51 |
| HALLMARK_ALLOGRAFT_REJECTION | 1.055 | 1.150 | -0.095 | 4.0e-49 | 2.6e-48 |


## Table S4a. Gene composition of the three a-priori nutrition–inflammation axes

16 genes per axis, entered as a-priori definition (not data-driven selection). All 48 symbols were detected in the expression matrix. 'detected_in_expression_matrix = N' would indicate a symbol dropped before scoring; there are none in this dataset.

| axis | order_in_set | gene_symbol | detected_in_expression_matrix |
|---|---|---|---|
| nutrition_liver | 1 | ALB | Y |
| nutrition_liver | 2 | TF | Y |
| nutrition_liver | 3 | RBP4 | Y |
| nutrition_liver | 4 | APOA1 | Y |
| nutrition_liver | 5 | APOA2 | Y |
| nutrition_liver | 6 | APOB | Y |
| nutrition_liver | 7 | GC | Y |
| nutrition_liver | 8 | TTR | Y |
| nutrition_liver | 9 | HPX | Y |
| nutrition_liver | 10 | SERPINA1 | Y |
| nutrition_liver | 11 | AHSG | Y |
| nutrition_liver | 12 | ORM1 | Y |
| nutrition_liver | 13 | FGA | Y |
| nutrition_liver | 14 | FGG | Y |
| nutrition_liver | 15 | ITIH2 | Y |
| nutrition_liver | 16 | PLG | Y |
| neutrophil | 1 | ELANE | Y |
| neutrophil | 2 | MPO | Y |
| neutrophil | 3 | LTF | Y |
| neutrophil | 4 | CEACAM8 | Y |
| neutrophil | 5 | S100A8 | Y |
| neutrophil | 6 | S100A9 | Y |
| neutrophil | 7 | S100A12 | Y |
| neutrophil | 8 | MMP9 | Y |
| neutrophil | 9 | PADI4 | Y |
| neutrophil | 10 | CXCR2 | Y |
| neutrophil | 11 | FPR1 | Y |
| neutrophil | 12 | DEFA4 | Y |
| neutrophil | 13 | CAMP | Y |
| neutrophil | 14 | LCN2 | Y |
| neutrophil | 15 | ARG1 | Y |
| neutrophil | 16 | OLFM4 | Y |
| lymphocyte | 1 | CD3D | Y |
| lymphocyte | 2 | CD3E | Y |
| lymphocyte | 3 | CD3G | Y |
| lymphocyte | 4 | CD8A | Y |
| lymphocyte | 5 | CD8B | Y |
| lymphocyte | 6 | LCK | Y |
| lymphocyte | 7 | GZMB | Y |
| lymphocyte | 8 | GZMK | Y |
| lymphocyte | 9 | PRF1 | Y |
| lymphocyte | 10 | IL7R | Y |
| lymphocyte | 11 | CCL5 | Y |
| lymphocyte | 12 | CXCL10 | Y |
| lymphocyte | 13 | CXCR6 | Y |
| lymphocyte | 14 | CD27 | Y |
| lymphocyte | 15 | CXCL13 | Y |
| lymphocyte | 16 | HLA-DRA | Y |


## Table S4b. MSigDB immune features used as pathway-level activity proxies

13 features (10 GO:BP + 3 Hallmark). 'set_size_in_gmt' is the raw MSigDB entry size; 'set_size_in_matrix' is the size after intersecting with the expression matrix. Overlap columns give the number of genes shared with each axis: the nutrition_liver axis shares 0 genes with all 13 features, which is what makes its correlations non-definitional. Source release is MSigDB v2024.1.Hs, established by byte-level matching against the public MSigDB release archive: md5 of h.all.symbols.gmt (d3606c325b9862bf3ea7991c2c2bf7d7) and of c5.go.bp.symbols.gmt (fd535c44825f9be3ea7c5d31e3d45f1d) each match v2024.1.Hs and no other release between v7.5.1 and v2026.1.Hs.

| immune_feature | feature_label | collection | set_size_in_gmt | set_size_in_matrix | overlap_nutrition_liver | overlap_neutrophil | overlap_lymphocyte | shared_with_any_axis |
|---|---|---|---|---|---|---|---|---|
| GOBP_T_CELL_ACTIVATION | T cell activation | GO:BP | 579 | 566 | 0 | 1 | 10 | Y |
| GOBP_B_CELL_ACTIVATION | B cell activation | GO:BP | 278 | 269 | 0 | 0 | 2 | Y |
| GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION | Antigen processing and presentation | GO:BP | 118 | 113 | 0 | 0 | 2 | Y |
| GOBP_MACROPHAGE_ACTIVATION | Macrophage activation | GO:BP | 112 | 109 | 0 | 0 | 0 | N |
| GOBP_NATURAL_KILLER_CELL_ACTIVATION | Natural killer cell activation | GO:BP | 98 | 97 | 0 | 0 | 0 | N |
| GOBP_MAST_CELL_ACTIVATION | Mast cell activation | GO:BP | 72 | 70 | 0 | 1 | 0 | Y |
| GOBP_REGULATORY_T_CELL_DIFFERENTIATION | Regulatory T cell differentiation | GO:BP | 41 | 40 | 0 | 0 | 1 | Y |
| GOBP_MYELOID_DENDRITIC_CELL_ACTIVATION | Myeloid dendritic cell activation | GO:BP | 29 | 29 | 0 | 0 | 0 | N |
| GOBP_NEUTROPHIL_DEGRANULATION | Neutrophil degranulation | GO:BP | 14 | 14 | 0 | 0 | 0 | N |
| GOBP_MONOCYTE_ACTIVATION | Monocyte activation | GO:BP | 12 | 12 | 0 | 0 | 0 | N |
| HALLMARK_INFLAMMATORY_RESPONSE | Inflammatory response | Hallmark | 200 | 199 | 0 | 1 | 5 | Y |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | Interferon gamma response | Hallmark | 200 | 196 | 0 | 1 | 2 | Y |
| HALLMARK_ALLOGRAFT_REJECTION | Allograft rejection | Hallmark | 200 | 195 | 0 | 2 | 11 | Y |


## Table S4c. Immune checkpoint and effector transcript panel

9 immune checkpoint genes and 10 effector genes, all detected in the expression matrix. Correlations for these transcripts are reported in Table 3 and Figure 2b.

| panel | order_in_panel | gene_symbol | detected_in_expression_matrix |
|---|---|---|---|
| Immune checkpoint | 1 | PDCD1 | Y |
| Immune checkpoint | 2 | CD274 | Y |
| Immune checkpoint | 3 | PDCD1LG2 | Y |
| Immune checkpoint | 4 | CTLA4 | Y |
| Immune checkpoint | 5 | LAG3 | Y |
| Immune checkpoint | 6 | HAVCR2 | Y |
| Immune checkpoint | 7 | TIGIT | Y |
| Immune checkpoint | 8 | IDO1 | Y |
| Immune checkpoint | 9 | SIGLEC15 | Y |
| Immune effector | 1 | IFNG | Y |
| Immune effector | 2 | GZMA | Y |
| Immune effector | 3 | GZMB | Y |
| Immune effector | 4 | PRF1 | Y |
| Immune effector | 5 | CXCL9 | Y |
| Immune effector | 6 | CXCL10 | Y |
| Immune effector | 7 | CXCL11 | Y |
| Immune effector | 8 | HLA-DRA | Y |
| Immune effector | 9 | STAT1 | Y |
| Immune effector | 10 | CCR5 | Y |


## Table S4d. MSigDB pathways used in the S1 feature scheme

6 real MSigDB pathways (source release v2024.1.Hs) spanning 14–732 genes after intersecting with the expression matrix (1,325 unique gene-pathway memberships in total). This heterogeneity in set size is why pathway scores were standardised within the training set before LASSO-Cox.

| pathway | collection | pathway_label | set_size_in_gmt | set_size_in_matrix |
|---|---|---|---|---|
| GOBP_NEUTROPHIL_DEGRANULATION | GO:BP | Neutrophil degranulation | 14 | 14 |
| GOBP_ACUTE_PHASE_RESPONSE | GO:BP | Acute phase response | 47 | 47 |
| HALLMARK_COAGULATION | Hallmark | Coagulation | 138 | 138 |
| HALLMARK_ALLOGRAFT_REJECTION | Hallmark | Allograft rejection | 200 | 195 |
| HALLMARK_INFLAMMATORY_RESPONSE | Hallmark | Inflammatory response | 200 | 199 |
| GOBP_ADAPTIVE_IMMUNE_RESPONSE | GO:BP | Adaptive immune response | 756 | 732 |


---

## 表格使用说明

- 主表：Table 1（基线）、Table 2（轴 × 免疫特征）、Table 3（轴 × 检查点/效应基因）、Table 4（内部验证，阴性）、Table 5（营养轴分组差异，独立）。
- 附表：S1 终点敏感性、S2 特征方案敏感性、S3 综合分分组（含循环论证警示）。
- **重点叙事**：Table 5 是本方案唯一不含循环论证的关联结果——营养轴与全部免疫特征零基因重叠，却与髓系/炎症通路活性系统性正相关；Table 4 如实报告预后模型无判别力。
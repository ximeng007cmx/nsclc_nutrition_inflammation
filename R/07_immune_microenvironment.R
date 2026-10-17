# =====================================================================
# 07_immune_microenvironment.R — Path A 转型：营养-炎症特征 × 免疫微环境
# 目标：不做预后模型的高判别力门槛，改问「营养-炎症分子特征与肿瘤免疫
#       微环境/免疫检查点是否相关」——弱信号也能产出信息量。
# 数据：全部来自本地方案（表达矩阵 + ssGSEA 打分 + MSigDB 本地 gmt），零新增下载。
# 说明：免疫特征为「通路级活性代理」（GO BP + Hallmark 真实条目，已核实存在），
#       非严格细胞去卷积；正文须如此表述，不得写成 CIBERSORT 估计值。
# 输出：
#   results/immune_axis_correlation.tsv  三轴 × 免疫特征 Spearman
#   results/checkpoint_correlation.tsv   三轴 × 免疫检查点基因
#   results/composite_group_diff.tsv     综合风险分分组（high vs low）免疫特征差异 → ESM 4 / Table S3
#                                         （原名 immune_group_diff.tsv，2026-10-17 改名以免与营养轴分组混淆）
#   figures/_intermediate/fig5_immune_heatmap.png
#   （【2026-10-16 A6-0.6】该热图为过程图：其信息已由稿件 Figure 2 覆盖，
#     不进投稿包；输出移至 figures/_intermediate/，原图内中文标题一并去除）
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

d   <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc3 <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
ex  <- d$expr_log
axes <- names(d$genesets)

# ---------- 1. 免疫特征基因集（全部为 MSigDB 真实条目，已 grep 核实） ----------
read_gmt <- function(path, want) {
  ln <- readLines(path)
  out <- list()
  for (l in ln) {
    p <- strsplit(l, "\t")[[1]]
    if (!(p[1] %in% want)) next
    out[[p[1]]] <- p[-c(1, 2)]
  }
  miss <- setdiff(want, names(out))
  if (length(miss)) stop("gmt 未找到: ", paste(miss, collapse = ", "))
  out
}
HALL <- file.path(PATHA$dir, "data/genesets/h.all.symbols.gmt")
GOBP <- file.path(PATHA$dir, "data/genesets/c5.go.bp.symbols.gmt")

IMM_GO <- c("GOBP_T_CELL_ACTIVATION", "GOBP_B_CELL_ACTIVATION",
            "GOBP_MACROPHAGE_ACTIVATION", "GOBP_NATURAL_KILLER_CELL_ACTIVATION",
            "GOBP_MYELOID_DENDRITIC_CELL_ACTIVATION",
            "GOBP_REGULATORY_T_CELL_DIFFERENTIATION",
            "GOBP_MAST_CELL_ACTIVATION", "GOBP_MONOCYTE_ACTIVATION",
            "GOBP_NEUTROPHIL_DEGRANULATION",
            "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION")
IMM_HALL <- c("HALLMARK_INFLAMMATORY_RESPONSE",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE",
              "HALLMARK_ALLOGRAFT_REJECTION")
imm_sets <- c(read_gmt(GOBP, IMM_GO), read_gmt(HALL, IMM_HALL))
imm_sets <- lapply(imm_sets, function(g) intersect(g, rownames(ex)))
imm_sets <- imm_sets[vapply(imm_sets, length, 0L) >= 10]
message("[07] 免疫特征 ", length(imm_sets), " 条；基因数范围 ",
        min(vapply(imm_sets, length, 0L)), "–", max(vapply(imm_sets, length, 0L)))

# 免疫检查点 / 效应分子面板（真实基因名）
CHECKPOINT <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2",
                "TIGIT", "IDO1", "SIGLEC15")
EFFECTOR   <- c("IFNG", "GZMA", "GZMB", "PRF1", "CXCL9", "CXCL10", "CXCL11",
                "HLA-DRA", "STAT1", "CCR5")
CK <- intersect(CHECKPOINT, rownames(ex))
EF <- intersect(EFFECTOR, rownames(ex))
message("[07] 检查点基因命中 ", length(CK), "/", length(CHECKPOINT),
        "；效应分子命中 ", length(EF), "/", length(EFFECTOR))

IMM_CACHE <- file.path(PATHA$work, "pa_scores_immune.rds")
if (file.exists(IMM_CACHE)) {
  imm <- readRDS(IMM_CACHE); message("[07] 免疫打分读缓存")
} else {
  imm <- ssgsea_score(ex, imm_sets)
  saveRDS(imm, IMM_CACHE)
}
imm_names <- colnames(imm)

# ---------- 2. 相关分析 ----------
sp <- function(x, y) {
  ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
  c(rho = unname(ct$estimate), p = ct$p.value)
}
S <- as.data.frame(sc3[, axes]); rownames(S) <- sc3$sample
common <- intersect(rownames(S), rownames(imm))
S <- S[common, , drop = FALSE]
I <- as.data.frame(imm[common, , drop = FALSE])

# 2a 三轴 × 免疫特征（附基因重叠度：区分「内部一致性」与「独立关联」）
r1 <- rbindlist(lapply(axes, function(a) {
  rbindlist(lapply(imm_names, function(b) {
    v <- sp(S[[a]], I[[b]])
    ov <- length(intersect(d$genesets[[a]], imm_sets[[b]]))
    data.table(axis = a, immune_feature = b, rho = v["rho"], p = v["p"],
               overlap_genes = ov, n_set = length(imm_sets[[b]]),
               kind = ifelse(ov > 0, "homologous (shared genes)",
                             "independent (no shared genes)"))
  }))
}))
r1[, fdr := p.adjust(p, "BH")]
fwrite(r1, file.path(PATHA$results, "immune_axis_correlation.tsv"), sep = "\t")
message("\n[07] 按「是否与免疫特征共享基因」分层的相关强度：")
print(r1[, .(n = .N, rho_median = round(median(rho), 3),
             rho_max = round(max(rho), 3)), by = .(axis, kind)][order(axis, kind)])

# 2b 三轴 × 免疫检查点/效应基因
gene_expr <- t(ex[c(CK, EF), common, drop = FALSE])
r2 <- rbindlist(lapply(axes, function(a) {
  rbindlist(lapply(colnames(gene_expr), function(g) {
    v <- sp(S[[a]], gene_expr[, g])
    data.table(axis = a, gene = g, rho = v["rho"], p = v["p"])
  }))
}))
r2[, fdr := p.adjust(p, "BH")]
fwrite(r2, file.path(PATHA$results, "checkpoint_correlation.tsv"), sep = "\t")

# 2c 高/低「营养-炎症综合分」组间免疫特征差异
#    综合分 = 三轴 z 分数按 03 模型系数线性组合（模型不可用时用等权 z）
mdl_f <- file.path(PATHA$work, "pa_model.rds")
if (file.exists(mdl_f)) {
  mdl <- readRDS(mdl_f)
  z <- as.data.frame(scale(S[, mdl$axes], center = mdl$mu, scale = mdl$sd))
  rs <- risk_score(z, mdl$coef)
} else {
  z <- as.data.frame(scale(S)); rs <- rowMeans(z)
}
grp <- ifelse(rs > median(rs), "High", "Low")
r3 <- rbindlist(lapply(imm_names, function(b) {
  x <- I[[b]]; w <- suppressWarnings(stats::wilcox.test(x ~ grp))
  data.table(immune_feature = b,
             median_high = median(x[grp == "High"]),
             median_low = median(x[grp == "Low"]),
             diff = median(x[grp == "High"]) - median(x[grp == "Low"]),
             p = w$p.value)
}))
r3[, fdr := p.adjust(p, "BH")]
setorder(r3, -diff)
fwrite(r3, file.path(PATHA$results, "composite_group_diff.tsv"), sep = "\t")

message("\n[07] 三轴与免疫特征相关（|rho| 前 8）：")
print(head(r1[order(-abs(rho)), .(axis, immune_feature, rho = round(rho, 3),
                                  p = signif(p, 2), fdr = signif(fdr, 2))], 8))
message("\n[07] 高低组免疫特征差异（P 前 8）：")
print(head(r3[order(p), .(immune_feature, diff = round(diff, 3),
                          p = signif(p, 2), fdr = signif(fdr, 2))], 8))

# ---------- 3. 热图 ----------
suppressPackageStartupMessages({ library(ggplot2); library(pheatmap) })
mat <- dcast(r1, immune_feature ~ axis, value.var = "rho")
m <- as.matrix(mat[, -1]); rownames(m) <- mat$immune_feature
m <- m[order(m[, 1]), , drop = FALSE]
INT <- file.path(PATHA$figures, "_intermediate")
dir.create(INT, showWarnings = FALSE, recursive = TRUE)
pheatmap(m,
         filename = file.path(INT, "fig5_immune_heatmap.png"),
         width = 5, height = 6, border_color = NA,
         color = grDevices::colorRampPalette(c("#2C6E9B", "white", "#C0392B"))(50))

message("[07] 完成 -> results/immune_axis_correlation.tsv 等 3 表 + fig5")

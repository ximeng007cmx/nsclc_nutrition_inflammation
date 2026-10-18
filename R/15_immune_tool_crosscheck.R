# =====================================================================
# 15_immune_tool_crosscheck.R — 独立免疫估计工具复核（评估第 2 项）
#
# 动机：主分析里 13 个免疫特征全部来自 MSigDB（GO:BP + Hallmark）经 ssGSEA 打分，
# 属同一工具族；审稿人可质疑「关联是 ssGSEA 的产物」。本脚本用**标记基因法**的
# MCP-counter（Becht et al., Genome Biol 2016）在**同一表达矩阵**上独立估计
# 10 个免疫/基质群体丰度，检验：
#   Q1  三条轴与各群体的关联方向是否与 ssGSEA 结论一致（尤其 nutrition 轴的髓系）
#   Q2  该方向是否在控制肿瘤纯度后仍成立（与 R/12 的偏相关同法）
#   Q3  两套工具本身是否互相印证（MCP 群体 × 13 个 ssGSEA 特征的中位相关）
#
# 依赖：MCPcounter 1.2.0 官方源码（work/MCPcounter-master/Source/R/MCPcounter.R，
#       因本机 R CMD INSTALL 段错误改为 source 调用）；签名文件
#       work/MCPcounter-master/Signatures/genes.txt，md5 944af31812c3f2aa61526074604e1c77
#
# 只读：work/pa_data.rds、work/pa_scores.rds、work/pa_scores_immune.rds、
#       results/estimate_scores_matched.tsv（不重算 ssGSEA）
#
# 输出：
#   results/immune_tool_crosscheck.tsv          轴 × MCP 群体（raw + 纯度偏相关）
#   results/immune_tool_crosscheck_summary.tsv  逐轴方向一致性汇总
#   results/immune_tool_concordance.tsv         MCP 群体 × 13 个 ssGSEA 特征（工具互证）
#   logs/run_15.log
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

log_con <- file(file.path(PATHA$logs, "run_15.log"), open = "wt")
sink(log_con, split = TRUE)
on.exit({ sink(); close(log_con) }, add = TRUE)

cat("== 15_immune_tool_crosscheck.R ==\n")

d   <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc  <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
imm <- readRDS(file.path(PATHA$work, "pa_scores_immune.rds"))
est <- fread(file.path(PATHA$results, "estimate_scores_matched.tsv"))

ex <- d$expr_log
axes <- names(d$genesets)
cat(sprintf("表达矩阵: %d 基因 × %d 样本；轴 %d 条；ssGSEA 免疫特征 %d 条\n",
            nrow(ex), ncol(ex), length(axes), ncol(imm)))

# ---- 1. MCP-counter（官方源码 + 本地签名）----
mcp_dir <- file.path(PATHA$work, "MCPcounter-master")
source(file.path(mcp_dir, "Source", "R", "MCPcounter.R"))
mcp_genes <- read.table(file.path(mcp_dir, "Signatures", "genes.txt"),
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                        check.names = FALSE, colClasses = "character")
mcp <- MCPcounter.estimate(ex, featuresType = "HUGO_symbols", genes = mcp_genes)
cat(sprintf("MCP-counter 群体 %d 个 × 样本 %d\n", nrow(mcp), ncol(mcp)))
print(rownames(mcp))

M <- as.data.frame(t(mcp), check.names = FALSE)
M$sample <- rownames(M)
ML <- c("Monocytic lineage", "Neutrophils", "Myeloid dendritic cells")

# ---- 2. 轴 × MCP 群体（Spearman + 纯度偏相关，与 R/12 同法）----
psp <- function(x, y, z) {
  rxy <- stats::cor(x, y, method = "spearman")
  rxz <- stats::cor(x, z, method = "spearman")
  ryz <- stats::cor(y, z, method = "spearman")
  den <- sqrt((1 - rxz^2) * (1 - ryz^2))
  r <- (rxy - rxz * ryz) / den
  r[abs(den) < 1e-12] <- NA_real_
  df <- length(x) - 3
  t <- r * sqrt(df / (1 - r^2))
  list(r = r, p = 2 * stats::pt(-abs(t), df))
}

S  <- sc[match(M$sample, sc$sample), axes, drop = FALSE]
sub <- est[matched == TRUE]
keep <- sub$barcode15 %in% M$sample
sub <- sub[keep]

rows <- rbindlist(lapply(axes, function(a) {
  x <- S[[a]]
  rbindlist(lapply(setdiff(colnames(M), "sample"), function(p) {
    y <- M[[p]]
    ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
    # 纯度偏相关子集
    idx <- match(sub$barcode15, M$sample)
    pp <- psp(x[idx], y[idx], sub$purity)
    data.table(axis = a, population = p, n = length(x),
               rho = unname(ct$estimate), p = ct$p.value,
               rho_adj_purity = pp$r, p_adj_purity = pp$p,
               myeloid = p %in% ML)
  }))
}))
rows[, `:=`(fdr = p.adjust(p, "BH"), fdr_adj_purity = p.adjust(p_adj_purity, "BH")),
     by = axis]
fwrite(rows, file.path(PATHA$results, "immune_tool_crosscheck.tsv"), sep = "\t")

# ---- 3. 工具互证：MCP 群体 × 13 个 ssGSEA 特征 ----
cm <- as.data.frame(imm[rownames(imm) %in% M$sample, , drop = FALSE])
M2 <- M[match(rownames(cm), M$sample), setdiff(colnames(M), "sample"), drop = FALSE]
conc <- rbindlist(lapply(colnames(M2), function(p) {
  rbindlist(lapply(colnames(cm), function(f) {
    ct <- suppressWarnings(stats::cor.test(M2[[p]], cm[[f]], method = "spearman"))
    data.table(population = p, ssgsea_feature = f,
               rho = unname(ct$estimate), p = ct$p.value)
  }))
}))
fwrite(conc, file.path(PATHA$results, "immune_tool_concordance.tsv"), sep = "\t")

# ---- 4. 汇总 ----
summ <- rows[, .(n_populations = .N,
                 n_positive = sum(rho > 0),
                 n_sig_positive = sum(rho > 0 & fdr < 0.05),
                 n_sig_negative = sum(rho < 0 & fdr < 0.05),
                 rho_median = median(rho),
                 myeloid_rho_median = median(rho[myeloid]),
                 myeloid_n_positive = sum(rho[myeloid] > 0),
                 myeloid_n_sig_positive = sum(rho[myeloid] > 0 & fdr[myeloid] < 0.05),
                 myeloid_rho_adj_median = median(rho_adj_purity[myeloid]),
                 myeloid_n_sig_adj = sum(rho_adj_purity[myeloid] > 0 &
                                         fdr_adj_purity[myeloid] < 0.05)),
             by = axis]
fwrite(summ, file.path(PATHA$results, "immune_tool_crosscheck_summary.tsv"), sep = "\t")

cat("\n-- 轴 × MCP-counter（合并 1,011 例）--\n")
print(rows[, .(axis, population, rho = round(rho, 3), fdr = signif(fdr, 2),
               rho_adj_purity = round(rho_adj_purity, 3))])
cat("\n-- 汇总 --\n")
print(summ)
cat("\n-- 工具互证（MCP 群体 × ssGSEA 特征，中位 rho）--\n")
print(conc[, .(median_rho = round(median(rho), 3),
               n_positive = sum(rho > 0), n = .N), by = population][order(-median_rho)])
cat("\n[15] 完成。\n")

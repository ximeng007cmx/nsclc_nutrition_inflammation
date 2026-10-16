# =====================================================================
# tools/audit_methods_facts.R — A3 Methods 写作前的只读事实审计
# 目的：把 Methods 中要写的每一个计数/参数，逐条从源产物实测取回，
#       禁止凭印象填写。不修改任何数据，只打印。
# 运行：Rscript tools/audit_methods_facts.R
# =====================================================================
suppressPackageStartupMessages(library(data.table))
setwd("D:/WorkBuddy_Y02/out_ali_directions/pathA")

d <- readRDS("work/pa_data.rds")

cat("---- 1. 表达矩阵 ----\n")
cat(sprintf("genes = %d ; samples = %d\n", nrow(d$expr_log), ncol(d$expr_log)))
cat(sprintf("min(expr_log) = %.4f ; log2(0.001) = %.4f\n",
            min(d$expr_log), log2(0.001)))
st <- sort(unique(substr(colnames(d$expr_log), 16, 17)))
cat(sprintf("sample type codes present = %s\n", paste(st, collapse = ",")))

cat("\n---- 2. 生存表 ----\n")
cat(sprintf("LUAD n=%d events=%d\n", nrow(d$surv_luad), sum(d$surv_luad$status == 1)))
cat(sprintf("LUSC n=%d events=%d\n", nrow(d$surv_lusc), sum(d$surv_lusc$status == 1)))
cat(sprintf("LUAD age NA=%d gender NA=%d stage_late NA=%d\n",
            sum(is.na(d$surv_luad$age)), sum(is.na(d$surv_luad$gender)),
            sum(is.na(d$surv_luad$stage_late))))
cat(sprintf("LUSC age NA=%d gender NA=%d stage_late NA=%d\n",
            sum(is.na(d$surv_lusc$age)), sum(is.na(d$surv_lusc$gender)),
            sum(is.na(d$surv_lusc$stage_late))))
cat("LUAD stage_late 分布：\n"); print(table(d$surv_luad$stage_late, useNA = "ifany"))
cat("LUSC stage_late 分布：\n"); print(table(d$surv_lusc$stage_late, useNA = "ifany"))

cat("\n---- 3. 三轴基因集 ----\n")
cat("geneset sizes: ",
    paste(sprintf("%s=%d", names(d$genesets), vapply(d$genesets, length, 0L)),
          collapse = ", "), "\n")
cat("missing genes : ",
    paste(sprintf("%s=%d", names(d$genesets_missing),
                  vapply(d$genesets_missing, length, 0L)), collapse = ", "), "\n")

cat("\n---- 4. 免疫特征（13 条）规模 ----\n")
immr <- readRDS("work/pa_scores_immune.rds")
cat(sprintf("immune features = %d ; samples = %d\n", ncol(immr), nrow(immr)))
cat("feature names:\n"); print(colnames(immr))

cat("\n---- 5. 相关性分析样本数与重叠分层 ----\n")
r1 <- fread("results/immune_axis_correlation.tsv")
cat(sprintf("rows = %d (3 axes x 13 features)\n", nrow(r1)))
print(r1[, .(n = .N, zero_overlap = sum(overlap_genes == 0),
             min_ov = min(overlap_genes), max_ov = max(overlap_genes)),
          by = axis])

cat("\n---- 6. 检查点/效应基因 ----\n")
r2 <- fread("results/checkpoint_correlation.tsv")
cat(sprintf("rows = %d ; unique genes = %d\n", nrow(r2), uniqueN(r2$gene)))
print(r2[, .(n = .N, n_sig = sum(p < 0.05), rho_min = round(min(rho), 3),
             rho_max = round(max(rho), 3)), by = axis])

cat("\n---- 7. Table 5 分组（营养轴中位数二分）----\n")
t5 <- fread("results/table5_nutrition_axis_raw.tsv")
cat(sprintf("features = %d ; n_high = %s ; n_low = %s\n",
            nrow(t5), unique(t5$n_high), unique(t5$n_low)))
cat(sprintf("all diff > 0 : %s ; all FDR < 0.05 : %s\n",
            all(t5$diff > 0), all(t5$fdr < 0.05)))
cat(sprintf("min P = %s ; n P underflowed to 0 = %d\n",
            format(min(t5$p), scientific = TRUE), sum(t5$p == 0)))

cat("\n---- 8. 内部验证 ----\n")
i1 <- fread("results/internal_split_summary.tsv")
print(i1)
rr <- fread("results/internal_repeated_cv.tsv")
cat(sprintf("repeats = %d ; median C_test = %s ; n distinct selected = %d\n",
            nrow(rr), round(median(rr$c_test), 4), uniqueN(rr$selected)))
print(rr[, .(med = round(median(c_test), 4),
             lo  = round(quantile(c_test, .025), 4),
             hi  = round(quantile(c_test, .975), 4),
             mean_C = round(mean(c_test), 4),
             med_opt = round(median(c_train - c_test), 4),
             mean_opt = round(mean(c_train - c_test), 4))])

cat("\n---- 9. 软硬件版本 ----\n")
cat(R.version.string, "\n")
for (p in c("data.table", "glmnet", "survival", "GSVA", "pheatmap", "ggplot2")) {
  cat(sprintf("  %-12s %s\n", p, as.character(packageVersion(p))))
}

cat("\n审计完成（只读）。\n")

# =====================================================================
# verify_table4.R —— 独立复核 Table 4（营养轴分组），并持久化未舍入原始值
# 2026-10-17 换号：原 Table 5 → Table 4（按正文首次引用顺序重排）；文件名一并改为 verify_table4.R
# ---------------------------------------------------------------------
# 目的：
#   1) 独立于 09_tables.R 重算一次 Table 5，确认结果可复现；
#   2) 把未舍入的原始值导出为 results/table5_nutrition_axis_raw.tsv，
#      供正文写作逐数字回查（Table 4 的 CSV 已四舍五入到 3 位小数）。
# 注意：
#   - results/composite_group_diff.tsv 是【综合风险分】分组的产物（对应 Table S3 / ESM 4），
#     其 diff 全为负属定义性结果，**不是** Table 5 的来源。勿混用。
#   - 本脚本只读 .rds、只写一个新 tsv，不修改任何既有产物。
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")
suppressPackageStartupMessages(library(data.table))

sc  <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
imm <- readRDS(file.path(PATHA$work, "pa_scores_immune.rds"))

S1 <- data.table(sample = sc$sample, nutrition_liver = sc$nutrition_liver)
IM <- as.data.table(imm)
IM[, sample := rownames(imm)]
mm <- merge(S1, IM, by = "sample")
mm[, grp := ifelse(nutrition_liver > median(nutrition_liver), "High", "Low")]

feat_cols <- setdiff(names(mm), c("sample", "nutrition_liver", "grp"))

res <- rbindlist(lapply(feat_cols, function(b) {
  x <- mm[[b]]; g <- mm$grp
  w <- suppressWarnings(stats::wilcox.test(x ~ g))
  data.table(
    immune_feature = b,
    n_high = sum(g == "High"),
    n_low  = sum(g == "Low"),
    median_high = median(x[g == "High"]),
    median_low  = median(x[g == "Low"]),
    diff        = median(x[g == "High"]) - median(x[g == "Low"]),
    p           = w$p.value
  )
}))
res[, fdr := p.adjust(p, "BH")]
setorder(res, -diff)

# 同时写两个文件名：R/10 与 manuscript 引用 table5_...，outline 记录 table4_...，
# 历史上二者内容相同（2026-10-17 数据修复后一并重生成）
OUT <- file.path(PATHA$results, "table4_nutrition_axis_raw.tsv")
OUT2 <- file.path(PATHA$results, "table5_nutrition_axis_raw.tsv")
fwrite(res, OUT, sep = "\t")
fwrite(res, OUT2, sep = "\t")

message(sprintf("[verify_table4] 样本数 n=%d（High %d / Low %d）",
                nrow(mm), sum(mm$grp == "High"), sum(mm$grp == "Low")))
message(sprintf("[verify_table4] 特征数=%d；diff>0 的 %d 个；FDR<0.05 的 %d 个",
                nrow(res), sum(res$diff > 0), sum(res$fdr < 0.05)))
message(sprintf("[verify_table4] 已写出未舍入原始值 -> %s", OUT))

# =====================================================================
# 02_score.R — 对 LUAD/LUSC 全部原发肿瘤样本做 ssGSEA 打分（3 条营养-炎症轴）
# 输出 work/pa_scores.rds：data.frame（样本 × 3 轴分数）
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

d <- readRDS(file.path(PATHA$work, "pa_data.rds"))
message(sprintf("[02] ssGSEA：%d 基因 × %d 样本，%d 个基因集",
                nrow(d$expr_log), ncol(d$expr_log), length(d$genesets)))
sc <- ssgsea_score(d$expr_log, d$genesets)
sc$sample <- rownames(sc)
sc <- sc[, c("sample", names(d$genesets)), drop = FALSE]
message("[02] 打分完成：")
print(summary(sc[, names(d$genesets)]))

saveRDS(sc, file.path(PATHA$work, "pa_scores.rds"))
message("[02] 完成 -> work/pa_scores.rds")

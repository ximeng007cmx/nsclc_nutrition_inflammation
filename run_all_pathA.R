# =====================================================================
# run_all_pathA.R — 路径 A 总控入口
# 用法（命令行）：
#   Rscript run_all_pathA.R            # 全流程
#   Rscript run_all_pathA.R 99         # 仅自测
#   Rscript run_all_pathA.R 01 02 03 04  # 选择性运行
#
# 【完整重建顺序】01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11
#   注意：11（Table S4）必须在生成表格索引之前跑完，因为索引会读取
#   results/tables/TableS4{a,b,c,d}_*.csv。索引生成：
#     python tools/make_article_tables.py
# =====================================================================
args <- commandArgs(trailingOnly = TRUE)
steps <- if (length(args) == 0) c("01", "02", "03", "04") else args

setwd("D:/WorkBuddy_Y02/out_ali_directions/pathA")
for (s in steps) {
  cat("\n========== 步骤", s, "==========\n")
  t0 <- Sys.time()
  source(file.path("R", paste0(s, switch(s,
    "01" = "_load.R",
    "02" = "_score.R",
    "03" = "_lasso_cox.R",
    "04" = "_validate.R",
    "05" = "_endpoint_sensitivity.R",
    "06" = "_geneset_sensitivity.R",
    "07" = "_immune_microenvironment.R",
    "08" = "_internal_validation.R",
    "09" = "_tables.R",
    "10" = "_figures_pathA.R",
    "11" = "_tableS4.R",
    "99" = "_selftest.R"))), echo = FALSE)
  cat(sprintf("---- 步骤 %s 完成，用时 %.1f 秒 ----\n",
              s, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# =====================================================================
# run_all_pathA.R — 路径 A 总控入口
# 用法（命令行）：
#   Rscript run_all_pathA.R            # 全流程
#   Rscript run_all_pathA.R 99         # 仅自测
#   Rscript run_all_pathA.R 01 02 03 04  # 选择性运行
#
# 【完整重建顺序】01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11 → 12 → 13
#   注意：11（Table S4）必须在生成表格索引之前跑完，因为索引会读取
#   results/tables/TableS4{a,b,c,d}_*.csv。索引生成：
#     python tools/make_article_tables.py
#   注意：12（随机面板零模型 + ESTIMATE 纯度校正）与 13（Fig. 6/7）依赖
#   data/external/estimate/ 下的三个 ESTIMATE 文件（见 DATA_SOURCES.md §8），
#   文件缺失时自动跳过并给出提示（其余步骤不受影响）。12 运行约 20 分钟。
# =====================================================================
args <- commandArgs(trailingOnly = TRUE)
steps <- if (length(args) == 0) c("01", "02", "03", "04") else args

setwd("D:/WorkBuddy_Y02/out_ali_directions/pathA")
for (s in steps) {
  cat("\n========== 步骤", s, "==========\n")
  t0 <- Sys.time()
  if (s %in% c("12", "13") &&
      !file.exists(file.path("data", "external", "estimate", "ESTIMATE_scores.txt"))) {
    cat("---- 步骤", s, "跳过：缺少 data/external/estimate/ESTIMATE_scores.txt",
        "（下载方法见 DATA_SOURCES.md §8）----\n")
    next
  }
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
    "12" = "_random_panel_estimate.R",
    "13" = "_figures_random_estimate.R",
    "99" = "_selftest.R"))), echo = FALSE)
  cat(sprintf("---- 步骤 %s 完成，用时 %.1f 秒 ----\n",
              s, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

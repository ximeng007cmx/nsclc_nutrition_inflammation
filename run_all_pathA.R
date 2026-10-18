# =====================================================================
# run_all_pathA.R — 路径 A 总控入口
# 用法（命令行）：
#   Rscript run_all_pathA.R            # 全流程
#   Rscript run_all_pathA.R 99         # 仅自测
#   Rscript run_all_pathA.R 01 02 03 04  # 选择性运行
#
# 【完整重建顺序】01 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11 → 12 → 13 → 14 → 15 → 16
#   注意：11（Table S4）必须在生成表格索引之前跑完，因为索引会读取
#   results/tables/TableS4{a,b,c,d}_*.csv。索引生成：
#     python tools/make_article_tables.py
#   注意：12（随机面板零模型 + ESTIMATE 纯度校正）与 13（Fig. 6/7）依赖
#   data/external/estimate/ 下的三个 ESTIMATE 文件（见 DATA_SOURCES.md §8），
#   文件缺失时自动跳过并给出提示（其余步骤不受影响）。12 运行约 20 分钟。
#   14（纯度/免疫含量增量 Cox）与 15（MCP-counter 独立复核）只读已落盘对象，秒级完成；
#   15 依赖 work/MCPcounter-master/（官方源码，见 `manuscript/immune_tool_crosscheck_2026-10-18.md`）。
#   16（独立外部队列验证）只读 data/external/{gse68465,gse31210}/ 下的
#   expr_symbol.tsv.gz 与 clinical.tsv（由 tools/fetch_external_tables.py →
#   tools/build_external_matrix.py → tools/build_external_clinical.py 生成，
#   原始数据来自 EBI BioStudies 的 ArrayExpress 镜像，见
#   `manuscript/external_validation_2026-10-18.md`）；缺目录时自动跳过。
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
  if (s == "16") {
    miss <- Filter(function(p) !file.exists(p), c(
      file.path("data", "external", "gse68465", "expr_symbol.tsv.gz"),
      file.path("data", "external", "gse68465", "clinical.tsv"),
      file.path("data", "external", "gse31210", "expr_symbol.tsv.gz"),
      file.path("data", "external", "gse31210", "clinical.tsv")))
    if (length(miss) > 0) {
      cat("---- 步骤 16 跳过：缺少", paste(miss, collapse = ", "),
          "（构建方法见 manuscript/external_validation_2026-10-18.md）----\n")
      next
    }
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
    "14" = "_incremental_cox.R",
    "15" = "_immune_tool_crosscheck.R",
    "16" = "_external_validation.R",
    "99" = "_selftest.R"))), echo = FALSE)
  cat(sprintf("---- 步骤 %s 完成，用时 %.1f 秒 ----\n",
              s, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

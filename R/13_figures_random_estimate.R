# R/13_figures_random_estimate.R
#
# 目的：把 R/12 的两项稳健性分析产物升格为正文主图（Fig. 6 / Fig. 7）。
# 输入（全部为 R/12 已落盘结果，不重算 ssGSEA / Cox / 零面板）：
#   - results/random_panel_draws.tsv.gz（K=1000 × 3 轴 × 2 零模型的 mean|rho| 抽样分布）
#   - results/estimate_partial_correlation.tsv（39 对 raw vs 纯度偏相关）
# 输出（A6-0 主图标准：宽 6.8 in ≤172.7 mm，PNG 600 dpi + cairo_pdf 矢量，图内无标题）：
#   - figures/fig6_random_null.png/.pdf
#   - figures/fig7_purity_adjust.png/.pdf
#
# 注：R/12 §5 仍会在 figures/_intermediate/ 生成 150 dpi 探索版（供人工判断），
#     投稿包一律以本脚本产物为准。

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
source("R/00_setup_pathA.R")
FIG <- PATHA$figures

axes <- c("lymphocyte", "neutrophil", "nutrition_liver")

# ---------- Fig. 6：随机面板零分布 vs 观测轴（mean |rho|） ----------
gz_path <- file.path(PATHA$results, "random_panel_draws.tsv.gz")
draws <- tryCatch(fread(gz_path), error = function(e)
  as.data.table(read.delim(gzfile(gz_path))))
draws <- draws[null_type != "observed", .(axis, null_type, value = rho_absmean, set = "null")]
obs <- draws[0][, set := "observed"]
# 观测值从零分布文件取（immune_rho_absmean，expr_matched 零模型行）
nsum <- fread(file.path(PATHA$results, "random_panel_null_summary.tsv"))
obs <- nsum[null_type == "expr_matched" & metric == "immune_rho_absmean",
            .(axis, null_type = "observed", value = obs, set = "observed")]
obs <- obs[, .(axis, null_type, value, set)]
obs$axis <- factor(obs$axis, levels = axes)
draws$axis <- factor(draws$axis, levels = axes)

pd <- rbind(draws, obs)
pd$null_type <- factor(pd$null_type,
                       levels = c("expr_matched", "unconstrained", "observed"),
                       labels = c("Expression-matched panels",
                                  "Unconstrained panels", "Observed axis"))
fig6 <- ggplot(pd, aes(x = value, fill = null_type)) +
  geom_histogram(bins = 40, alpha = .6, position = "identity",
                 colour = "grey30", linewidth = .1) +
  geom_vline(data = pd[set == "observed"], aes(xintercept = value),
             colour = "#C0392B", linewidth = .9) +
  facet_wrap(~axis, scales = "free_y",
             labeller = labeller(axis = c(lymphocyte = "Lymphocyte axis",
                                          neutrophil = "Neutrophil axis",
                                          nutrition_liver = "Nutrition axis"))) +
  labs(x = expression("Mean |Spearman " * rho * "| across 13 immune features"),
       y = "Random panels",
       fill = NULL) +
  scale_fill_manual(values = c("Expression-matched panels" = "#5B8DB8",
                               "Unconstrained panels" = "#A9C6DE",
                               "Observed axis" = "#C0392B")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", strip.background = element_rect(fill = "grey92"))
ggsave(file.path(FIG, "fig6_random_null.png"), fig6,
       width = 6.8, height = 3.8, dpi = 600)
ggsave(file.path(FIG, "fig6_random_null.pdf"), fig6,
       width = 6.8, height = 3.8, device = cairo_pdf)

# ---------- Fig. 7：纯度校正前后（39 对散点，按基因共享分层） ----------
part <- fread(file.path(PATHA$results, "estimate_partial_correlation.tsv"))
part$axis <- factor(part$axis, levels = axes)
part$kind <- factor(part$kind,
                    levels = c("independent (no shared genes)",
                               "homologous (shared genes)"),
                    labels = c("No shared genes", "Shared genes"))
fig7 <- ggplot(part, aes(x = rho_raw, y = rho_adj_purity, colour = kind, shape = kind)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 1.7, alpha = .85) +
  facet_wrap(~axis,
             labeller = labeller(axis = c(lymphocyte = "Lymphocyte axis",
                                          neutrophil = "Neutrophil axis",
                                          nutrition_liver = "Nutrition axis"))) +
  labs(x = expression("Spearman " * rho * " (raw)"),
       y = expression("Partial " * rho * " (purity-adjusted)"),
       colour = NULL, shape = NULL) +
  scale_colour_manual(values = c("No shared genes" = "#1F6FB4",
                                 "Shared genes" = "#C0392B")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", strip.background = element_rect(fill = "grey92"))
ggsave(file.path(FIG, "fig7_purity_adjust.png"), fig7,
       width = 6.8, height = 3.4, dpi = 600)
ggsave(file.path(FIG, "fig7_purity_adjust.pdf"), fig7,
       width = 6.8, height = 3.4, device = cairo_pdf)

message("\n[13] Fig. 6 / Fig. 7 已写出 -> figures/fig6_random_null.{png,pdf}, fig7_purity_adjust.{png,pdf}")
message("[13] 观测 mean|rho|："); print(obs)

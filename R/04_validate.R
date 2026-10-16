# =====================================================================
# 04_validate.R — 三级验证：LUAD 内部（bootstrap C-index / KM / 多变量 Cox）
#                 + LUSC 外部泛化
# 输出：
#   results/validation_summary.tsv
#   figures/supplementary/FigureS1_km_luad.png + .pdf
#   figures/supplementary/FigureS2_km_lusc.png + .pdf
#   figures/supplementary/FigureS3_lasso_lambda.png + .pdf
#   figures/supplementary/FigureS4_forest_luad.png + .pdf
# 【2026-10-16 A6-0.1/0.3/0.4 修改】
#   - 附图由本脚本直接写入 figures/supplementary/ 并带最终文件名（此前写
#     figures/fig{1,2,3,4}_*.png 再手工拷贝，重跑不会更新附图，属可复现性断点，
#     已修复；旧顶层文件移至 figures/_intermediate/legacy_20261016/）
#   - 按 Discover Oncology 指南：图内不放标题（"Do not include titles or
#     captions within your illustrations"），log-rank P 等信息在图注中；
#     PNG 600 dpi，另出 cairo_pdf 矢量版（字体嵌入）
#   - 只改输出路径与样式，不改任何计算逻辑与随机数
# 注意：验证集高低风险组切点沿用训练集中位数（避免乐观偏倚）。
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

d   <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc  <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
mdl <- readRDS(file.path(PATHA$work, "pa_model.rds"))

zfun <- function(df) as.data.frame(scale(df[, mdl$axes, with = FALSE],
                                        center = mdl$mu, scale = mdl$sd))
prep <- function(sv) {
  dt <- merge(sv, sc, by = "sample")
  dt$risk_score <- risk_score(zfun(dt), mdl$coef)
  dt$risk_grp <- ifelse(dt$risk_score > mdl$train_median, "High", "Low")
  dt
}
tr <- prep(d$surv_luad)
va <- prep(d$surv_lusc)

summary_rows <- list()

# ---- LUAD 训练集 ----
safe_logrank <- function(time, status, grp) {
  if (length(unique(grp)) < 2) return(NA_real_)
  logrank(time, status, grp)
}
safe_km_dt <- function(dt, cohort) {
  if (length(unique(dt$risk_grp)) < 2) return(NULL)
  km_dt(dt, cohort)
}

ci_tr <- cindex_boot(tr$time, tr$status, tr$risk_score,
                     n = PARS$boot_n, seed = PARS$seed)
lr_tr <- safe_logrank(tr$time, tr$status, tr$risk_grp)
summary_rows$luad <- data.table(cohort = "TCGA-LUAD (train)", n = nrow(tr),
                                events = sum(tr$status == 1),
                                cindex = ci_tr["est"], ci_lo = ci_tr["lo"],
                                ci_hi = ci_tr["hi"], logrank_p = lr_tr)

# ---- LUSC 外部泛化 ----
ci_va <- cindex_boot(va$time, va$status, va$risk_score,
                     n = PARS$boot_n, seed = PARS$seed)
lr_va <- safe_logrank(va$time, va$status, va$risk_grp)
summary_rows$lusc <- data.table(cohort = "TCGA-LUSC (external)", n = nrow(va),
                                events = sum(va$status == 1),
                                cindex = ci_va["est"], ci_lo = ci_va["lo"],
                                ci_hi = ci_va["hi"], logrank_p = lr_va)

# ---- LUAD 多变量 Cox（风险分 + 年龄 + 性别 + 分期）----
mv <- cox_multivariate(
  tr[!is.na(age) & !is.na(gender) & !is.na(stage_late)],
  covars = c("age", "gender", "stage_late"))
summary_rows$mv <- mv$table
fwrite(mv$table, file.path(PATHA$results, "multivariate_cox_luad.tsv"), sep = "\t")

out <- rbind(summary_rows$luad, summary_rows$lusc)
fwrite(out, file.path(PATHA$results, "validation_summary.tsv"), sep = "\t")

# ---- 图 ----
suppressPackageStartupMessages(library(ggplot2))

km_dt <- function(dt, cohort) {
  f <- survival::survfit(survival::Surv(time, status) ~ risk_grp, data = dt)
  sf <- summary(f, times = sort(unique(dt$time)))
  g <- sub("^risk_grp=", "", as.character(sf$strata))
  data.table(time = sf$time, surv = sf$surv, grp = g, cohort = cohort)
}
draw_km <- function(dt, file_stem) {
  p <- ggplot(dt, aes(time, surv, colour = grp)) +
    geom_step(linewidth = 0.8) +
    scale_colour_manual(values = c(High = "#C0392B", Low = "#1F6FB4")) +
    labs(x = "Time (days)", y = "Overall survival", colour = "Risk") +
    theme_bw(base_size = 12)
  ggsave(paste0(file_stem, ".png"), p, width = 5.5, height = 4, dpi = 600)
  ggsave(paste0(file_stem, ".pdf"), p, width = 5.5, height = 4, device = cairo_pdf)
}
SUP <- file.path(PATHA$figures, "supplementary")
km_luad <- safe_km_dt(tr, "LUAD")
km_lusc <- safe_km_dt(va, "LUSC")
if (!is.null(km_luad))
  draw_km(km_luad, file.path(SUP, "FigureS1_km_luad"))
if (!is.null(km_lusc))
  draw_km(km_lusc, file.path(SUP, "FigureS2_km_lusc"))

# LASSO lambda 选择图
p3 <- ggplot(mdl$cv, aes(lambda, cvm)) +
  geom_errorbar(aes(ymin = cvm - cvsd, ymax = cvm + cvsd), alpha = 0.3) +
  geom_point(size = 0.8) +
  geom_vline(xintercept = mdl$lambda_min, linetype = 2, colour = "#C0392B") +
  scale_x_log10() +
  labs(x = expression(lambda), y = "Partial likelihood deviance") +
  theme_bw(base_size = 12)
ggsave(file.path(SUP, "FigureS3_lasso_lambda.png"), p3,
       width = 5.5, height = 4, dpi = 600)
ggsave(file.path(SUP, "FigureS3_lasso_lambda.pdf"), p3,
       width = 5.5, height = 4, device = cairo_pdf)

# 森林图（多变量 Cox）
tab <- mv$table
tab[, term := factor(term, levels = rev(term))]
p4 <- ggplot(tab, aes(HR, term)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_x_log10() +
  labs(x = "Hazard ratio (log scale)", y = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(SUP, "FigureS4_forest_luad.png"), p4,
       width = 5.5, height = 4, dpi = 600)
ggsave(file.path(SUP, "FigureS4_forest_luad.pdf"), p4,
       width = 5.5, height = 4, device = cairo_pdf)

message("[04] 三级验证完成：")
print(out)
print(mv$table)

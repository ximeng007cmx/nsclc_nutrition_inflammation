# =====================================================================
# 08_internal_validation.R — LUAD 内部 train/test 验证（替代跨组织学外部验证）
# 依据 06 的实测结论：跨组织学（LUAD↔LUSC）外部 C 天花板约 0.55，
#   且特征越多内部 C 越高、外部越低（过拟合）。故改为：
#   ① 简单随机 7:3 分割 LUAD（`sample.int()`，**未**按事件分层；
#      实测训练集事件率 38.3% vs 测试集 32.0%，分层本会给 36.4%）；
#   ② 特征选择与标准化**只在训练集**内完成（防泄漏）；
#   ③ 测试集评估；④ 重复 100 次随机分割报告分布（比单次分割稳健）。
# 输出：
#   results/internal_split_summary.tsv     单次分割（含训练/测试 C）
#   results/internal_repeated_cv.tsv       100 次重复分割的测试集 C 分布
#   figures/_intermediate/fig6_internal_km.png
#   （【2026-10-16 A6-0.6】该 KM 为过程图：同信息已由稿件 Figure 5c 复现，
#     不进投稿包；输出移至 figures/_intermediate/，原图内标题一并去除）
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

d   <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc3 <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
sv  <- copy(d$surv_luad)
sv  <- merge(sv, sc3, by = "sample")
axes <- names(d$genesets)
message(sprintf("[08] LUAD 可用样本 %d（事件 %d）", nrow(sv), sum(sv$status == 1)))

fit_once <- function(tr, te, p_univ = 0.10) {
  mu <- vapply(tr[, ..axes], mean, 0); sdv <- vapply(tr[, ..axes], stats::sd, 0)
  ztr <- as.data.frame(scale(tr[, ..axes], center = mu, scale = sdv))
  zte <- as.data.frame(scale(te[, ..axes], center = mu, scale = sdv))
  # 训练集内单变量 Cox 初筛
  pv <- vapply(axes, function(a) {
    summary(survival::coxph(survival::Surv(tr$time, tr$status) ~ ztr[[a]]))$coefficients[1, 5]
  }, 0)
  keep <- names(pv)[pv <= p_univ]
  if (!length(keep)) keep <- names(pv)[which.min(pv)]   # 至少保留一个
  if (length(keep) == 1L) {
    # 仅 1 个特征无法做 LASSO（glmnet 要求 >=2 列）→ 用单变量 Cox 系数
    f1 <- survival::coxph(survival::Surv(tr$time, tr$status) ~ ztr[[keep]])
    cf <- setNames(unname(stats::coef(f1)), keep)
    if (!is.finite(cf) || cf == 0) cf <- setNames(1e-8, keep)
  } else {
    fit <- suppressWarnings(fit_lasso_cox(ztr[, keep, drop = FALSE], tr$time, tr$status,
                                          nfolds = PARS$lasso_nfolds, seed = PARS$seed))
    cf <- fit$coef[fit$coef != 0]
    if (!length(cf)) cf <- setNames(1e-8, keep[1])
  }
  rtr <- risk_score(ztr[, names(cf), drop = FALSE], cf)
  rte <- risk_score(zte[, names(cf), drop = FALSE], cf)
  cut <- median(rtr)
  list(cf = cf, keep = keep, rtr = rtr, rte = rte, cut = cut,
       c_tr = cindex(tr$time, tr$status, rtr),
       c_te = if (stats::sd(rte) > 0) cindex(te$time, te$status, rte) else 0.5,
       p_te = safe_lr(te$time, te$status, rte, cut))
}
safe_lr <- function(time, status, rs, cut) {
  g <- ifelse(rs > cut, "High", "Low")
  if (length(unique(g)) < 2) return(NA_real_)
  logrank(time, status, g)
}

# ---------- 单次分割（简单随机 7:3，未按事件分层） ----------
set.seed(PARS$seed)
idx <- sample.int(nrow(sv), floor(0.7 * nrow(sv)))
tr <- sv[idx]; te <- sv[-idx]
message(sprintf("[08] 单次分割：训练 %d（事件 %d）/ 测试 %d（事件 %d）",
                nrow(tr), sum(tr$status == 1), nrow(te), sum(te$status == 1)))
r1 <- fit_once(tr, te)
out1 <- data.table(
  scheme = "LUAD internal 7:3",
  n_train = nrow(tr), n_test = nrow(te),
  selected = paste(names(r1$cf), collapse = "+"),
  coef = paste(round(r1$cf, 4), collapse = "+"),
  c_train = round(r1$c_tr, 4), c_test = round(r1$c_te, 4),
  p_test = signif(r1$p_te, 3))
fwrite(out1, file.path(PATHA$results, "internal_split_summary.tsv"), sep = "\t")
message("[08] 单次分割结果："); print(out1)

# ---------- 重复 100 次分割（每次均为简单随机 7:3，非分层） ----------
B <- 100
rr <- rbindlist(lapply(seq_len(B), function(b) {
  set.seed(PARS$seed + b)
  ix <- sample.int(nrow(sv), floor(0.7 * nrow(sv)))
  o <- suppressWarnings(fit_once(sv[ix], sv[-ix]))
  data.table(rep = b, c_train = o$c_tr, c_test = o$c_te, n_sel = length(o$cf))
}))
summ <- rr[, .(
  c_train_median = round(median(c_train), 4),
  c_test_median = round(median(c_test), 4),
  c_test_q025 = round(quantile(c_test, .025), 4),
  c_test_q975 = round(quantile(c_test, .975), 4),
  frac_sel_gt0 = round(mean(n_sel > 0), 3))]
rr[, gap := c_train - c_test]
fwrite(rr, file.path(PATHA$results, "internal_repeated_cv.tsv"), sep = "\t")
message(sprintf("[08] 重复 %d 次分割：训练 C 中位 %.3f / 测试 C 中位 %.3f [%.3f, %.3f]；平均过拟合差 %.3f",
                B, summ$c_train_median, summ$c_test_median,
                summ$c_test_q025, summ$c_test_q975, mean(rr$gap)))

# ---------- KM（测试集） ----------
suppressPackageStartupMessages(library(ggplot2))
te2 <- copy(te)
te2[, rs := r1$rte][, grp := ifelse(rs > r1$cut, "High", "Low")]
if (length(unique(te2$grp)) >= 2) {
  f <- survival::survfit(survival::Surv(time, status) ~ grp, data = te2)
  sf <- summary(f, times = sort(unique(te2$time)))
  kd <- data.table(time = sf$time, surv = sf$surv,
                   grp = sub("^grp=", "", as.character(sf$strata)))
  p <- ggplot(kd, aes(time, surv, colour = grp)) +
    geom_step(linewidth = 0.8) +
    scale_colour_manual(values = c(High = "#C0392B", Low = "#1F6FB4")) +
    labs(x = "Time (days)", y = "Overall survival", colour = "Risk") +
    theme_bw(base_size = 12)
  INT <- file.path(PATHA$figures, "_intermediate")
  dir.create(INT, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(INT, "fig6_internal_km.png"), p,
         width = 5.5, height = 4, dpi = 600)
}
message("[08] 完成 -> results/internal_*.tsv + figures/_intermediate/fig6_internal_km.png")

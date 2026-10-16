# =====================================================================
# tools/audit_screen_pvalues.R — 只读诊断：初筛 P 值分布
# 目的：核实「100 次重复中每次恰有 1 个轴通过 P<0.10」是否合理，
#       还是初筛环节存在实现问题。打印每轴的 P 与最小值。
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

d   <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc3 <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
sv  <- merge(copy(d$surv_luad), sc3, by = "sample")
axes <- names(d$genesets)

pv_of <- function(tr) {
  mu <- vapply(tr[, ..axes], mean, 0); sdv <- vapply(tr[, ..axes], stats::sd, 0)
  ztr <- as.data.frame(scale(tr[, ..axes], center = mu, scale = sdv))
  vapply(axes, function(a) {
    summary(survival::coxph(survival::Surv(tr$time, tr$status) ~ ztr[[a]]))$coefficients[1, 5]
  }, 0)
}

cat("\n===== 全样本（n=500）单变量 P =====\n")
print(round(pv_of(sv), 4))

cat("\n===== 100 次 70% 子样本：每轴 P 的分位数 =====\n")
M <- t(vapply(seq_len(100), function(b) {
  set.seed(PARS$seed + b)
  ix <- sample.int(nrow(sv), floor(0.7 * nrow(sv)))
  pv_of(sv[ix])
}, numeric(length(axes))))
colnames(M) <- axes
for (a in axes) {
  q <- quantile(M[, a], c(0, .25, .5, .75, 1))
  cat(sprintf("  %-16s min=%.4f q25=%.4f median=%.4f q75=%.4f max=%.4f  |  P<=0.10 占比=%.0f%%\n",
              a, q[1], q[2], q[3], q[4], q[5], 100 * mean(M[, a] <= 0.10)))
}
cat(sprintf("\n  每次「至少一个轴 P<=0.10」的占比 = %.0f%%\n",
            100 * mean(apply(M, 1, min) <= 0.10)))
cat(sprintf("  每次「通过初筛的轴数」分布：\n"))
print(table(apply(M, 1, function(r) sum(r <= 0.10))))

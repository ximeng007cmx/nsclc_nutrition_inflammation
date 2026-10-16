# =====================================================================
# tools/audit_internal_selected.R — 只读复核 R/08 内部验证的两个待证事实
#
# 待证事实（A2 正文中已写入，需证据支持）：
#   F1. 分割是否为「按事件分层」？
#   F2. 100 次重复中，被选中的特征名是否每次都相同？
#       （results/internal_repeated_cv.tsv 只存 n_sel，未存特征名）
#
# 方法：按 R/08 完全相同的随机数构造（set.seed(PARS$seed + b)；
#       sample.int(n, floor(0.7n))）与 fit_once 逻辑重算，
#       并把 c_train / c_test 与已存结果逐行比对以证明复现一致。
# 输出：results/_audit_internal_selected.tsv（新文件，不覆盖任何既有结果）
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

d   <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc3 <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
sv  <- copy(d$surv_luad)
sv  <- merge(sv, sc3, by = "sample")
axes <- names(d$genesets)
message(sprintf("[audit] LUAD 可用样本 %d（事件 %d）", nrow(sv), sum(sv$status == 1)))

safe_lr <- function(time, status, rs, cut) {
  g <- ifelse(rs > cut, "High", "Low")
  if (length(unique(g)) < 2) return(NA_real_)
  logrank(time, status, g)
}

fit_once <- function(tr, te, p_univ = 0.10) {
  mu <- vapply(tr[, ..axes], mean, 0); sdv <- vapply(tr[, ..axes], stats::sd, 0)
  ztr <- as.data.frame(scale(tr[, ..axes], center = mu, scale = sdv))
  zte <- as.data.frame(scale(te[, ..axes], center = mu, scale = sdv))
  pv <- vapply(axes, function(a) {
    summary(survival::coxph(survival::Surv(tr$time, tr$status) ~ ztr[[a]]))$coefficients[1, 5]
  }, 0)
  n_pass <- sum(pv <= p_univ)                 # 真正通过初筛的轴数（兜底前）
  pass_names <- paste(names(pv)[pv <= p_univ], collapse = "+")
  keep <- names(pv)[pv <= p_univ]
  fellback <- !length(keep)                   # 无轴通过 → 是否走了「最小 P」兜底
  if (fellback) keep <- names(pv)[which.min(pv)]
  if (length(keep) == 1L) {
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
  list(cf = cf, keep = keep, cut = cut,
       n_pass = n_pass, pass_names = pass_names, fellback = fellback,
       pv = pv,
       c_tr = cindex(tr$time, tr$status, rtr),
       c_te = if (stats::sd(rte) > 0) cindex(te$time, te$status, rte) else 0.5,
       p_te = safe_lr(te$time, te$status, rte, cut))
}

# ---------- F1：分割是否分层？ ----------
set.seed(PARS$seed)
idx <- sample.int(nrow(sv), floor(0.7 * nrow(sv)))
tr <- sv[idx]; te <- sv[-idx]
n_tr <- nrow(tr); e_tr <- sum(tr$status == 1)
n_te <- nrow(te); e_te <- sum(te$status == 1)
cat("\n===== F1 单次分割的事件构成 =====\n")
cat(sprintf("训练 %d（事件 %d，%.1f%%） / 测试 %d（事件 %d，%.1f%%） / 全体事件率 %.1f%%\n",
            n_tr, e_tr, 100 * e_tr / n_tr, n_te, e_te, 100 * e_te / n_te,
            100 * sum(sv$status == 1) / nrow(sv)))
cat(sprintf("分层分割的期望测试事件数 = %.1f；实测 = %d\n",
            0.3 * sum(sv$status == 1), e_te))
cat("→ 结论：", if (abs(e_te - 0.3 * sum(sv$status == 1)) < 2) "与分层一致" else
    "与分层不一致，实际为简单随机分割（未分层）", "\n")

o1 <- fit_once(tr, te)
cat(sprintf("\n单次分割：训练集内单变量 P = %s\n",
            paste(sprintf("%s=%.4f", names(o1$pv), o1$pv), collapse = ", ")))
cat(sprintf("  通过 P<0.10 的轴数 = %d（%s）；是否走最小 P 兜底 = %s\n",
            o1$n_pass, if (nzchar(o1$pass_names)) o1$pass_names else "无", o1$fellback))
cat(sprintf("  入模特征 = %s；系数 = %s\n",
            paste(names(o1$cf), collapse = "+"),
            paste(round(o1$cf, 4), collapse = "+")))

# ---------- F2：100 次重复中选中的特征名 ----------
B <- 100
rr <- rbindlist(lapply(seq_len(B), function(b) {
  set.seed(PARS$seed + b)
  ix <- sample.int(nrow(sv), floor(0.7 * nrow(sv)))
  o <- suppressWarnings(fit_once(sv[ix], sv[-ix]))
  data.table(rep = b, c_train = o$c_tr, c_test = o$c_te,
             n_sel = length(o$cf),
             selected = paste(names(o$cf), collapse = "+"),
             n_pass_screen = o$n_pass,
             fellback = o$fellback)
}))
fwrite(rr, file.path(PATHA$results, "_audit_internal_selected.tsv"), sep = "\t")

cat("\n===== F2 100 次重复中入模的特征 =====\n")
print(rr[, .N, by = selected][order(-N)])
cat("\n入模特征名不同的次数 =", uniqueN(rr$selected), "\n")
cat("通过初筛的轴数分布：\n"); print(table(rr$n_pass_screen))
cat(sprintf("无轴通过初筛、走最小 P 兜底的次数 = %d / %d\n", sum(rr$fellback), B))

# ---------- 与已存结果逐行比对，证明复现一致 ----------
old <- fread(file.path(PATHA$results, "internal_repeated_cv.tsv"))
setorder(old, rep); setorder(rr, rep)
dtr <- max(abs(old$c_train - rr$c_train))
dte <- max(abs(old$c_test - rr$c_test))
cat("\n===== 复现一致性 =====\n")
cat(sprintf("max|c_train 差| = %.3g ; max|c_test 差| = %.3g ; n_sel 完全相同 = %s\n",
            dtr, dte, identical(old$n_sel, rr$n_sel)))
cat(if (dtr < 1e-12 && dte < 1e-12) "→ 复现精确一致，本次审计的特征名可代表 R/08 的结果\n" else
    "→ 复现不一致，审计结论不可用于替代 R/08 结果\n")

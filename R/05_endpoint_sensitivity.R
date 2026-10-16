# =====================================================================
# 05_endpoint_sensitivity.R — 终点敏感性：OS / DSS / PFI 三终点平行比较
# 目的：回答「切 PFI 是否真的能救信号？」——用同一管道、同一基因集重拟合。
# 输入：work/pa_data.rds（含 OS/DSS/PFI 三套 time/status）、work/pa_scores.rds
# 输出：results/endpoint_sensitivity.tsv（训练/外部 × 3 终点的 C-index 与 log-rank）
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

d  <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
axes <- names(d$genesets)

make_sv <- function(meta_proj, endpoint) {
  # 复用 01 已合并的临床（surv_* 表里已带三终点列），仅替换 time/status
  sv <- copy(meta_proj)
  sv[, time := get(paste0(endpoint, "_time"))]
  sv[, status := get(endpoint)]
  sv <- sv[!is.na(time) & !is.na(status) & time > 0]
  sv$time <- as.numeric(sv$time); sv$status <- as.numeric(sv$status)
  sv
}

rows <- list()
for (ep in c("OS", "DSS", "PFI")) {
  for (co in c("LUAD", "LUSC")) {
    sv <- make_sv(if (co == "LUAD") d$surv_luad else d$surv_lusc, ep)
    dt <- merge(sv, sc, by = "sample")
    if (nrow(dt) < 50 || sum(dt$status == 1) < 20) {
      rows[[paste(ep, co)]] <- data.table(endpoint = ep, cohort = co,
                                          n = nrow(dt), events = sum(dt$status == 1),
                                          cindex = NA_real_, logrank_p = NA_real_,
                                          n_axis_sel = 0L, note = "insufficient samples or events")
      next
    }
    mu <- vapply(dt[, ..axes], mean, 0); sdv <- vapply(dt[, ..axes], stats::sd, 0)
    z <- as.data.frame(scale(dt[, ..axes], center = mu, scale = sdv))
    fit <- suppressWarnings(
      fit_lasso_cox(z, dt$time, dt$status, nfolds = PARS$lasso_nfolds,
                    seed = PARS$seed))
    rs <- risk_score(z, fit$coef)
    if (all(fit$coef == 0) || stats::sd(rs) == 0) {
      ci <- 0.5; lr <- NA_real_
    } else {
      ci <- cindex(dt$time, dt$status, rs)
      lr <- logrank(dt$time, dt$status, ifelse(rs > median(rs), "High", "Low"))
    }
    # 单变量最强轴的 P（用于判断是否有任何轴单独显著）
    pmin_axis <- min(vapply(axes, function(a) {
      summary(survival::coxph(survival::Surv(dt$time, dt$status) ~ z[[a]]))$coefficients[1, 5]
    }, 0))
    rows[[paste(ep, co)]] <- data.table(
      endpoint = ep, cohort = co, n = nrow(dt), events = sum(dt$status == 1),
      cindex = ci, logrank_p = lr,
      n_axis_sel = sum(fit$coef != 0),
      note = sprintf("strongest single axis P=%.3g; s_used=%s", pmin_axis,
                     ifelse(is.character(fit$s_used), fit$s_used,
                            sprintf("%.4g", fit$s_used))))
  }
}
out <- rbindlist(rows)
out[, cindex := round(cindex, 4)][, logrank_p := round(logrank_p, 4)]
fwrite(out, file.path(PATHA$results, "endpoint_sensitivity.tsv"), sep = "\t")
message("[05] 终点敏感性结果：")
print(out)

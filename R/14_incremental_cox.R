# =====================================================================
# 14_incremental_cox.R — 纯度/免疫含量的多变量 Cox 增量检验
#
# 目的（对应 2026-10-18 评估第 1 项）：把「轴分/复合分测的是纯度与免疫含量，
# 不是独立预后生物学」从相关层面（R/12 的偏相关）推进到回归层面：
#    Q1  risk_score 在校正 purity / ImmuneScore 后是否仍显著？
#    Q2  purity / ImmuneScore 加入临床模型是否有增量（LRT）？
#    Q3  反过来，risk_score 在 clinical + purity + ImmuneScore 之上是否有增量？
#    Q4  各模型 C-index 与增量 ΔC-index
#
# 队列：TCGA-LUAD 生存队列（与 04_validate.R 的 multivariate_cox_luad 完全同一集合）
# 协变量：age + gender + stage_late（临床基座），purity（ESTIMATE，每 0.1），
#         ImmuneScore（每 SD）
# 只读已有落盘对象（pa_data / pa_scores / pa_model / estimate_scores_matched），
# 不重算打分与模型。
#
# 输出：
#   results/incremental_cox_terms.tsv    逐模型逐项 HR
#   results/incremental_cox_summary.tsv  逐模型 n/events/C-index/ΔC
#   results/incremental_cox_tests.tsv    似然比检验（增量是否成立）
#   logs/run_14.log
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

log_con <- file(file.path(PATHA$logs, "run_14.log"), open = "wt")
sink(log_con, split = TRUE)
on.exit({ sink(); close(log_con) }, add = TRUE)

cat("== 14_incremental_cox.R ==\n")

d   <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc  <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
mdl <- readRDS(file.path(PATHA$work, "pa_model.rds"))
est <- fread(file.path(PATHA$results, "estimate_scores_matched.tsv"))

# ---- 1. 组装 LUAD 分析集（与 04 的 prep() 完全一致，再并入 ESTIMATE）----
zfun <- function(df) as.data.frame(scale(df[, mdl$axes, with = FALSE],
                                        center = mdl$mu, scale = mdl$sd))
prep <- function(sv) {
  dt <- merge(sv, sc, by = "sample")
  dt$risk_score <- risk_score(zfun(dt), mdl$coef)
  dt
}
dt <- as.data.table(prep(d$surv_luad))
dt <- merge(dt, est[, .(sample = barcode15, purity, Immune, Stromal)],
            by = "sample", all.x = TRUE)
dt <- dt[!is.na(age) & !is.na(gender) & !is.na(stage_late) &
         !is.na(time) & !is.na(status) & !is.na(purity)]
dt[, `:=`(purity_10 = purity * 10,
          immune_z  = as.numeric(scale(Immune)),
          stromal_z = as.numeric(scale(Stromal)))]
# 单轴对照（lymphocyte 为全队列唯一单因素显著轴），按 SD 标准化
dt[, lymph_z := as.numeric(scale(lymphocyte))]

cat(sprintf("分析集: n = %d, events = %d, median OS = %.0f d\n",
            nrow(dt), sum(dt$status == 1), median(dt$time)))

# ---- 2. 模型族 ----
S <- function(rhs) {
  f <- stats::as.formula(paste("survival::Surv(time, status) ~", rhs))
  environment(f) <- globalenv()
  survival::coxph(f, data = dt, ties = "efron")
}

CLIN <- "age + gender + stage_late"
m0 <- S(CLIN)                                   # 临床基线（不含风险分）
m1 <- S(paste("risk_score +", CLIN))
m2 <- S(paste("risk_score +", CLIN, "+ purity_10"))
m3 <- S(paste("risk_score +", CLIN, "+ immune_z"))
m4 <- S(paste("risk_score +", CLIN, "+ purity_10 + immune_z"))
m5 <- S(paste(CLIN, "+ purity_10 + immune_z"))   # 无 risk_score

# 单轴版（lymphocyte，每 SD）
m1b <- S(paste("lymph_z +", CLIN))
m4b <- S(paste("lymph_z +", CLIN, "+ purity_10 + immune_z"))
# 注：轴对照块的参照模型与 M5 完全相同（同一公式），故不重复定义

models <- list(
  "M0 clinical only"                               = m0,
  "M1 clinical + risk_score"                       = m1,
  "M2 M1 + purity"                                 = m2,
  "M3 M1 + ImmuneScore"                            = m3,
  "M4 M1 + purity + ImmuneScore"                   = m4,
  "M5 clinical + purity + ImmuneScore (ref)"       = m5,

  "M1b clinical + lymphocyte axis (per SD)"        = m1b,
  "M4b M1b + purity + ImmuneScore"                 = m4b)

# ---- 3. 逐项 HR ----
term_rows <- rbindlist(lapply(names(models), function(nm) {
  fit <- models[[nm]]
  s <- summary(fit)
  data.table(model = nm, term = rownames(s$coefficients),
             HR = s$conf.int[, "exp(coef)"],
             lo = s$conf.int[, "lower .95"],
             hi = s$conf.int[, "upper .95"],
             p = s$coefficients[, "Pr(>|z|)"])
}))
fwrite(term_rows, file.path(PATHA$results, "incremental_cox_terms.tsv"), sep = "\t")

# ---- 4. 逐模型 C-index（bootstrap CI）----
ci_rows <- rbindlist(lapply(names(models), function(nm) {
  fit <- models[[nm]]
  lp <- as.numeric(predict(fit, type = "lp"))
  ci <- cindex_boot(dt$time, dt$status, lp, n = PARS$boot_n, seed = PARS$seed)
  data.table(model = nm, n = fit$n, events = fit$nevent,
             cindex = ci["est"], lo = ci["lo"], hi = ci["hi"],
             loglik = as.numeric(fit$loglik[2]),
             n_par = length(coef(fit)))
}))
c_ref <- ci_rows[grepl("^M5 ", model), cindex]
ci_rows[, delta_c := cindex - c_ref]
fwrite(ci_rows, file.path(PATHA$results, "incremental_cox_summary.tsv"), sep = "\t")

# ---- 5. 似然比检验 ----
lrt <- function(a, b, label, question) {
  an <- anova(a, b, test = "LRT")
  data.table(comparison = label, question = question,
             lr_chisq = an$Chisq[2], df = an$Df[2], p = an$`Pr(>|Chi|)`[2])
}
tests <- rbind(
  lrt(m0, m1,  "M0 -> M1 (+risk_score)",                 "risk_score adds to clinical covariates alone?"),
  lrt(m0, m1b, "M0 -> M1b (+lymphocyte axis)",          "lymphocyte axis adds to clinical covariates alone?"),
  lrt(m1, m2,  "M1 -> M2 (+purity)",                     "purity adds to clinical + risk_score?"),
  lrt(m1, m3,  "M1 -> M3 (+ImmuneScore)",                "ImmuneScore adds to clinical + risk_score?"),
  lrt(m1, m4,  "M1 -> M4 (+purity +ImmuneScore)",        "both add to clinical + risk_score?"),
  lrt(m0, m5,  "M0 -> M5 (+purity +ImmuneScore)",        "content adds to clinical covariates alone?"),
  lrt(m5, m4,  "M5 -> M4 (+risk_score)",                 "risk_score adds beyond purity + ImmuneScore?"),
  lrt(m5, m4b, "M5 -> M4b (+lymphocyte axis)",           "lymphocyte axis adds beyond purity + ImmuneScore?"))
tests[, p_fmt := formatC(p, format = "e", digits = 2)]
fwrite(tests, file.path(PATHA$results, "incremental_cox_tests.tsv"), sep = "\t")

# ---- 5b. 共线性：risk_score 与 ImmuneScore（解释衰减机制用）----
coll <- data.table(
  pair = c("risk_score_vs_ImmuneScore", "lymphocyte_vs_ImmuneScore",
           "risk_score_vs_purity", "ImmuneScore_vs_purity"),
  spearman_rho = c(stats::cor(dt$risk_score, dt$Immune, method = "spearman"),
                   stats::cor(dt$lymphocyte, dt$Immune, method = "spearman"),
                   stats::cor(dt$risk_score, dt$purity, method = "spearman"),
                   stats::cor(dt$Immune, dt$purity, method = "spearman")),
  n = nrow(dt))
fwrite(coll, file.path(PATHA$results, "incremental_cox_collinearity.tsv"), sep = "\t")

cat("\n-- 共线性 --\n")
print(coll[, .(pair, spearman_rho = round(spearman_rho, 4))])

cat("\n-- risk_score / lymphocyte HR 随校正变化 --\n")
key <- term_rows[term %in% c("risk_score", "lymph_z"),
                 .(model, term, HR = round(HR, 3), lo = round(lo, 3),
                   hi = round(hi, 3), p = signif(p, 3))]
print(key)

cat("\n-- C-index --\n")
print(ci_rows[, .(model, n, events, cindex = round(cindex, 3),
                  lo = round(lo, 3), hi = round(hi, 3),
                  delta_c = round(delta_c, 3))])

cat("\n-- LRT --\n")
print(tests[, .(comparison, lr_chisq = round(lr_chisq, 2), df, p = signif(p, 3))])

cat("\n[14] 完成。\n")

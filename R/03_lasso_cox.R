# =====================================================================
# 03_lasso_cox.R — 训练集（TCGA-LUAD）LASSO-Cox 建模
# 输出：
#   work/pa_model.rds   模型对象（系数、lambda、训练中位风险分）
#   results/model_coef.tsv
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

d  <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc <- readRDS(file.path(PATHA$work, "pa_scores.rds"))

tr <- merge(d$surv_luad, sc, by = "sample")
message(sprintf("[03] LUAD 训练表：%d 例（事件 %d）", nrow(tr), sum(tr$status == 1)))

# 三轴分数标准化（用训练集均值/SD，验证集沿用同一变换）
axes <- names(d$genesets)
mu <- vapply(tr[, ..axes], mean, 0)
sd <- vapply(tr[, ..axes], stats::sd, 0)
zfun <- function(df) as.data.frame(scale(df[, ..axes], center = mu, scale = sd))
tr_z <- zfun(tr)

message("[03] 单变量 Cox 诊断（三轴，LUAD OS）：")
uni <- do.call(rbind, lapply(axes, function(a) {
  s <- summary(survival::coxph(survival::Surv(tr$time, tr$status) ~ tr_z[[a]]))
  data.table(axis = a, HR = s$conf.int[1], lo = s$conf.int[3],
             hi = s$conf.int[4], p = s$coefficients[1, 5])
}))
print(uni)
fwrite(uni, file.path(PATHA$results, "axis_univariate_cox_luad.tsv"), sep = "\t")

message("[03] LASSO-Cox（10 折 CV）…")
fit <- fit_lasso_cox(tr_z, tr$time, tr$status,
                     nfolds = PARS$lasso_nfolds, alpha = PARS$lasso_alpha,
                     seed = PARS$seed)
if (fit$s_used != "lambda.min") {
  message(sprintf("[03][注意] lambda.min 全零，回退 lambda = %.4g（结果文件已标注）",
                  fit$s_used))
}
print(round(fit$coef, 4))

tr$risk_score <- risk_score(tr_z, fit$coef)
train_median <- median(tr$risk_score)

res <- data.table(axis = names(fit$coef), coef = fit$coef,
                  lambda_min = fit$lambda_min, s_used = fit$s_used)
fwrite(res, file.path(PATHA$results, "model_coef.tsv"), sep = "\t")

saveRDS(list(coef = fit$coef, lambda_min = fit$lambda_min,
             s_used = fit$s_used,
             mu = mu, sd = sd, train_median = train_median,
             axes = axes,
             cv = data.table(lambda = fit$cvfit$lambda,
                             cvm = fit$cvfit$cvm,
                             cvsd = fit$cvfit$cvsd)),
        file.path(PATHA$work, "pa_model.rds"))
message(sprintf("[03] 完成。风险分中位数（训练集切点）= %.4f", train_median))

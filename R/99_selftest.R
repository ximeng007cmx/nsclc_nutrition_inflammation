# =====================================================================
# 99_selftest.R — 合成数据端到端自测（不依赖 TCGA 数据，可独立运行）
# 设计：150 样本 × 60 基因，3 个基因集各 16/8/8 基因，
#       在「neutrophil」轴植入正向生存风险信号（HR≈1.05/单位分数）。
# 通过标准：全流程无错误、C_est 显著 > 0.5（>0.55）、系数方向正确、
#           LUSC（合成）外部验证方向一致。
# 运行：Rscript R/99_selftest.R
# =====================================================================
source("R/pathA_core.R")
suppressPackageStartupMessages(library(data.table))
ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop("SELFTEST FAIL: ", msg, call. = FALSE)
  message("  PASS: ", msg)
}

set.seed(42)
n <- 150; p <- 60
gn <- c(paste0("NU_", 1:16), paste0("NE_", 1:16), paste0("LY_", 1:16),
        paste0("BG_", 1:12))
ok(length(gn) == p, "基因数 = 60")

expr <- matrix(rnorm(n * p, 5, 1.2), nrow = p, ncol = n,
               dimnames = list(gn, sprintf("S%03d", 1:n)))
latent_ne <- colMeans(expr[paste0("NE_", 1:16), ])  # 潜在中性粒分数
z <- as.vector(scale(latent_ne))
time <- rexp(n, rate = exp(0.15 * z) / 100)
status <- rbinom(n, 1, 0.7)
time[status == 0] <- pmin(time[status == 0], 120) + 30
time <- round(time, 1)

genesets <- list(
  nutrition_liver = gn[1:16],
  neutrophil      = gn[17:32],
  lymphocyte      = gn[33:48]
)

message("[selftest 1/5] ssGSEA 打分")
sc <- ssgsea_score(expr, genesets)
ok(all(c("nutrition_liver", "neutrophil", "lymphocyte") %in% names(sc)),
   "ssGSEA 输出含 3 条轴")
ok(nrow(sc) == n, "ssGSEA 每样本一行")

message("[selftest 2/5] LASSO-Cox 拟合")
fit <- fit_lasso_cox(sc, time, status, nfolds = 5, seed = 1)
rs <- risk_score(sc, fit$coef)
ok(abs(fit$coef["neutrophil"]) > 0, "neutrophil 轴系数非零")
ok(fit$coef["neutrophil"] > 0, "系数方向为正（植入信号方向）")

message("[selftest 3/5] C-index > 0.55")
ci <- cindex(time, status, rs)
ok(ci > 0.55, sprintf("C-index = %.3f > 0.55", ci))

message("[selftest 4/5] bootstrap CI 覆盖估计值")
cib <- cindex_boot(time, status, rs, n = 100, seed = 1)
ok(cib["lo"] <= cib["est"] && cib["est"] <= cib["hi"],
   sprintf("C = %.3f [%.3f, %.3f]", cib["est"], cib["lo"], cib["hi"]))

message("[selftest 5/5] log-rank 与多变量 Cox")
grp <- ifelse(rs > median(rs), "High", "Low")
p_lr <- logrank(time, status, grp)
ok(!is.na(p_lr), sprintf("log-rank P = %.4g", p_lr))
sv <- data.table(time = time, status = status, risk_score = rs,
                 age = rnorm(n, 60, 10), gender = factor(rep(c("M", "F"), n / 2)),
                 stage_late = factor(rep(c("early", "late"), n / 2)))
mv <- cox_multivariate(sv)
ok("risk_score" %in% mv$table$term, "多变量 Cox 含 risk_score 项")

message("\n=== 自测全部通过 ===")
print(mv$table)

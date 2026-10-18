# =====================================================================
# 16_external_validation.R — 独立外部队列验证「阴性」（评估第 3 项）
#
# 队列（均为 NCBI 不可达时经 EBI BioStudies 获取的 ArrayExpress 镜像）：
#   GSE68465（E-GEOD-68465，Director's Challenge，LUAD，HG-U133A）442 例 / 236 事件
#   GSE31210（E-GEOD-31210，日本 LUAD IA–II，HG-U133 Plus 2）226 例 / 35 事件
#
# 设计：轴定义、ssGSEA 参数与风险分系数全部**冻结自 TCGA 训练**（mdl$coef, mdl$mu, mdl$sd），
#       不重新拟合、不重新选特征；评价指标与 R/04 的外部验证完全一致
#       （Harrell C-index + bootstrap CI，切点沿用训练集中位数）。
# 输入：
#   data/external/<study>/expr_symbol.tsv.gz（log2(x+1)，symbol × 样本，tools/build_external_matrix.py）
#   data/external/<study>/clinical.tsv（tools/build_external_clinical.py）
# 输出：
#   results/external_validation.tsv
#   results/external_axis_coverage.tsv
#   logs/run_16.log
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

log_con <- file(file.path(PATHA$logs, "run_16.log"), open = "wt")
sink(log_con, split = TRUE)
on.exit({ sink(); close(log_con) }, add = TRUE)

cat("== 16_external_validation.R ==\n")

d   <- readRDS(file.path(PATHA$work, "pa_data.rds"))
mdl <- readRDS(file.path(PATHA$work, "pa_model.rds"))
gs  <- d$genesets

COHORTS <- list(
  list(id = "GSE68465", dir = "gse68465", platform = "HG-U133A (GPL96)", histology = "LUAD"),
  list(id = "GSE31210", dir = "gse31210", platform = "HG-U133 Plus 2 (GPL570)", histology = "LUAD"))

cov_rows <- list(); val_rows <- list()

for (co in COHORTS) {
  base <- file.path(PATHA$dir, "data", "external", co$dir)
  ex <- as.matrix(fread(file.path(base, "expr_symbol.tsv.gz")), rownames = 1)
  cl <- fread(file.path(base, "clinical.tsv"), skip = 1)
  cl <- cl[!is.na(os_time_days) & os_time_days > 0]
  keep <- intersect(cl$sample, colnames(ex))
  cl <- cl[match(keep, sample)]
  ex <- ex[, keep, drop = FALSE]
  cat(sprintf("\n%s: %d 基因 × %d 样本（%d 事件）\n", co$id, nrow(ex), ncol(ex), sum(cl$os_event == 1)))

  # 轴基因覆盖
  for (a in names(gs)) {
    g <- intersect(gs[[a]], rownames(ex))
    cov_rows[[length(cov_rows) + 1]] <- data.table(
      cohort = co$id, axis = a, n_axis_genes = length(gs[[a]]),
      n_present = length(g), pct_present = round(100 * length(g) / length(gs[[a]]), 1))
  }

  # 与原分析同参数打分（仅在缺失基因时按可用基因打分，覆盖 < 100% 已在表中列明）
  sc <- ssgsea_score(ex, gs)
  sc <- as.data.frame(sc[keep, , drop = FALSE])

  # 冻结系数的风险分：用 TCGA 的 mu/sd 标准化，再乘 TCGA 系数
  z <- sweep(sweep(as.matrix(sc[, mdl$axes, drop = FALSE]), 2, mdl$mu, "-"), 2, mdl$sd, "/")
  rs <- as.vector(z[, names(mdl$coef), drop = FALSE] %*% mdl$coef)

  ci <- cindex_boot(cl$os_time_days, cl$os_event, rs, n = PARS$boot_n, seed = PARS$seed)
  # 切点：训练集中位数在 TCGA 内为 0（分数以 0 为中心），跨平台不可迁移，
  # 故 log-rank 用**队列内中位数**切分，并在输出中标注 cut 列；主指标为阈值无关的 C-index。
  grp <- ifelse(rs > stats::median(rs), "High", "Low")
  lr <- if (length(unique(grp)) > 1) logrank(cl$os_time_days, cl$os_event, grp) else NA_real_

  # 队列内单轴 Cox（每 SD），用于说明"换队列后单轴也不显著"
  ax_rows <- rbindlist(lapply(mdl$axes, function(a) {
    x <- as.numeric(scale(sc[[a]]))
    f <- survival::coxph(survival::Surv(cl$os_time_days, cl$os_event) ~ x)
    s <- summary(f)
    data.table(cohort = co$id, axis = a,
               HR_per_SD = s$conf.int[1, "exp(coef)"],
               lo = s$conf.int[1, "lower .95"], hi = s$conf.int[1, "upper .95"],
               p = s$coefficients[1, "Pr(>|z|)"])
  }))
  ax_rows[, fdr := p.adjust(p, "BH")]

  val_rows[[co$id]] <- rbind(
    data.table(cohort = co$id, platform = co$platform, histology = co$histology,
               analysis = "composite (TCGA coefficients, frozen)", n = ncol(ex),
               events = sum(cl$os_event == 1),
               cindex = ci["est"], lo = ci["lo"], hi = ci["hi"], logrank_p = lr,
               cut = "cohort median",
               HR = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_, p = NA_real_, fdr = NA_real_),
    ax_rows[, .(cohort, platform = co$platform, histology = co$histology,
                analysis = paste("single axis:", axis), n = ncol(ex),
                events = sum(cl$os_event == 1), cindex = NA_real_, lo = NA_real_,
                hi = NA_real_, logrank_p = NA_real_, cut = NA_character_,
                HR = HR_per_SD, HR_lo = lo, HR_hi = hi, p = p, fdr = fdr)])
}

val <- rbindlist(val_rows, fill = TRUE)
cov <- rbindlist(cov_rows)
fwrite(val, file.path(PATHA$results, "external_validation.tsv"), sep = "\t")
fwrite(cov, file.path(PATHA$results, "external_axis_coverage.tsv"), sep = "\t")

cat("\n-- 轴基因覆盖 --\n"); print(cov)
cat("\n-- 外部验证 --\n")
print(val[, .(cohort, analysis, n, events,
              cindex = ifelse(is.na(cindex), NA, round(cindex, 3)),
              lo = ifelse(is.na(lo), NA, round(lo, 3)),
              hi = ifelse(is.na(hi), NA, round(hi, 3)),
              logrank_p = signif(logrank_p, 3),
              HR = ifelse(is.na(HR), NA, round(HR, 3)),
              p = signif(p, 3))])
cat("\n[16] 完成。\n")

# =====================================================================
# 09_tables.R — 生成文章级表格（Table 1–5 + Table S1–S2）
# 数据来源：全部为已跑通的中间产物（work/*.rds、results/*.tsv），不新增分析。
# 输出：results/tables/*.csv（每个表独立，便于直接投稿或转 Word）
# 规范：P 值统一格式；比率保留 1 位小数；rho 保留 3 位小数。
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

TDIR <- file.path(PATHA$results, "tables")
dir.create(TDIR, showWarnings = FALSE, recursive = TRUE)

fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), "NA",
    ifelse(p < 1e-300, "<1e-300",
      ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
             formatC(p, format = "f", digits = digits))))
}
f3 <- function(x) formatC(x, format = "f", digits = 3)
f1 <- function(x) formatC(x, format = "f", digits = 1)

d   <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc  <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
axes <- names(d$genesets)

prep <- function(sv) as.data.table(merge(sv, sc, by = "sample"))
tr <- prep(d$surv_luad); va <- prep(d$surv_lusc)
message(sprintf("[09] 分析样本：LUAD %d（事件 %d）/ LUSC %d（事件 %d）",
                nrow(tr), sum(tr$status == 1), nrow(va), sum(va$status == 1)))

# ---------------- Table 1 队列基线特征 ----------------
iqr <- function(x, dig = 0) {
  q <- stats::quantile(x, c(.25, .5, .75), na.rm = TRUE)
  fmt <- paste0("%.", dig, "f (%.", dig, "f–%.", dig, "f)")
  sprintf(fmt, q[2], q[1], q[3])
}
st <- function(x) {
  x <- tolower(as.character(x))
  ifelse(grepl("^stage i{1,2}[abc]?$|^stage i{1,2}$", x), "Early (I–II)",
   ifelse(grepl("^stage ii{1,2}|^stage iv", x), "Late (III–IV)",
    ifelse(is.na(x) | x == "na" | x == "", "Unknown", "Other")))
}
baseline <- function(dt, label) {
  s <- table(st(dt$stage))
  cnt <- function(k) as.integer(ifelse(k %in% names(s), s[[k]], 0L))
  pct <- function(k) 100 * cnt(k) / nrow(dt)
  cens <- dt[status == 0]
  pctxt <- function(n, p) sprintf("%d (%.1f)", n, p)
  data.table(
    Characteristic = c(
      "Patients, n", "Deaths, n (%)", "Age, median (IQR)", "Male, n (%)",
      "Stage early (I–II), n (%)", "Stage late (III–IV), n (%)",
      "Stage other/unclassified, n (%)", "Stage unknown, n (%)",
      "OS time, median (IQR), days",
      "Follow-up of censored, median (IQR), days",
      "nutrition_liver score, median (IQR)",
      "neutrophil score, median (IQR)",
      "lymphocyte score, median (IQR)"),
    Value = c(
      as.character(nrow(dt)),
      pctxt(sum(dt$status == 1), 100 * mean(dt$status == 1)),
      iqr(dt$age),
      pctxt(sum(dt$gender == "MALE"), 100 * mean(dt$gender == "MALE")),
      pctxt(cnt("Early (I–II)"), pct("Early (I–II)")),
      pctxt(cnt("Late (III–IV)"), pct("Late (III–IV)")),
      pctxt(cnt("Other"), pct("Other")),
      pctxt(cnt("Unknown"), pct("Unknown")),
      iqr(dt$time),
      if (nrow(cens)) iqr(cens$time) else "NA",
      iqr(dt$nutrition_liver, 3), iqr(dt$neutrophil, 3), iqr(dt$lymphocyte, 3)))
}
t1 <- merge(baseline(tr, "LUAD"), baseline(va, "LUSC"),
            by = "Characteristic", suffixes = c("_LUAD", "_LUSC"))
t1 <- t1[, .(Characteristic, LUAD = Value_LUAD, LUSC = Value_LUSC)]
fwrite(t1, file.path(TDIR, "Table1_baseline.csv"))

# ---------------- Table 2 三轴 × 免疫特征 ----------------
r1 <- fread(file.path(PATHA$results, "immune_axis_correlation.tsv"))
setorder(r1, axis, -rho)          # 先按数值排序
r1[, type_en := ifelse(overlap_genes > 0, "Overlapping genes", "No overlapping genes")]
t2 <- r1[, .(axis, immune_feature, n_set, overlap_genes,
             rho = f3(rho), p = fmt_p(p), FDR = fmt_p(fdr), kind = type_en)]
setnames(t2, c("axis", "immune_feature", "n_set", "overlap_genes",
               "rho", "p", "FDR", "kind"),
         c("Nutrition–inflammation axis", "Immune feature (MSigDB)",
           "Set size", "Overlapping genes", "Spearman rho", "P", "FDR", "Type"))
fwrite(t2, file.path(TDIR, "Table2_axis_immune_correlation.csv"))

# ---------------- Table 3 三轴 × 免疫检查点/效应基因 ----------------
r2 <- fread(file.path(PATHA$results, "checkpoint_correlation.tsv"))
r2[, panel := ifelse(gene %in% c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3",
                                 "HAVCR2","TIGIT","IDO1","SIGLEC15"),
                     "Immune checkpoint", "Effector / TIS-like")]
t3r <- dcast(r2, gene + panel ~ axis, value.var = "rho")
pc  <- dcast(r2, gene ~ axis, value.var = "p")
t3 <- data.table(gene = t3r$gene, panel = t3r$panel)
for (a in axes) t3[[paste0(a, " (rho)")]] <- f3(t3r[[a]])
for (a in axes) t3[[paste0(a, " P")]] <- fmt_p(pc[[a]])
setorder(t3, panel, gene)
fwrite(t3, file.path(TDIR, "Table3_axis_checkpoint_correlation.csv"))

# ---------------- Table 5 内部验证 ----------------
# 2026-10-17 编号互换：内部验证在正文 §4 才首次引用，按「首次引用顺序」排为 Table 5
i1 <- fread(file.path(PATHA$results, "internal_split_summary.tsv"))
rr <- fread(file.path(PATHA$results, "internal_repeated_cv.tsv"))
t4 <- data.table(
  Analysis = c("Single 7:3 split — training set",
               "Single 7:3 split — held-out test set",
               sprintf("Repeated 100 splits — test set, median [2.5th–97.5th]"),
               "Mean optimism (training − test C)"),
  Value = c(f3(i1$c_train), sprintf("%s (log-rank P = %s)", f3(i1$c_test), fmt_p(i1$p_test)),
            sprintf("%s [%s–%s]", f3(median(rr$c_test)),
                    f3(quantile(rr$c_test, .025)), f3(quantile(rr$c_test, .975))),
            f3(mean(rr$c_train - rr$c_test))),
  Note = c("Features selected within training set only",
           "Two axes passed P<0.10 screening (lymphocyte P = 0.019, nutrition P = 0.053); both retained after LASSO shrinkage",
           "Same procedure re-run on 100 random splits",
           "Shows the degree of over-optimism in internally fitted C-index"))
fwrite(t4, file.path(TDIR, "Table5_internal_validation.csv"))

# ---------------- Table 4 按「营养轴」分组的免疫特征差异（无循环论证） ----------------
# 2026-10-17 编号互换：本节在正文 §3 首次引用，按「首次引用顺序」排为 Table 4
# 说明：不用综合风险分分组（其含 lymphocyte 轴负系数 → 差异必然为负，属同义反复），
#       改用 nutrition_liver 轴自身的中位数分组；该轴与全部免疫特征零基因重叠。
imm <- readRDS(file.path(PATHA$work, "pa_scores_immune.rds"))
S1 <- data.table(sample = sc$sample, nutrition_liver = sc$nutrition_liver)
IM <- as.data.table(imm); IM[, sample := rownames(imm)]
mm <- merge(S1, IM, by = "sample")
mm[, grp := ifelse(nutrition_liver > median(nutrition_liver), "High", "Low")]
feat_cols <- setdiff(names(mm), c("sample", "nutrition_liver", "grp"))
t5 <- rbindlist(lapply(feat_cols, function(b) {
  x <- mm[[b]]; g <- mm$grp
  w <- suppressWarnings(stats::wilcox.test(x ~ g))
  data.table(immune_feature = b,
             `Median (nutrition-high)` = median(x[g == "High"]),
             `Median (nutrition-low)` = median(x[g == "Low"]),
             Difference = median(x[g == "High"]) - median(x[g == "Low"]),
             P = w$p.value)
}))
t5[, FDR := p.adjust(P, "BH")]
setorder(t5, -Difference)
t5out <- t5[, .(immune_feature,
                `Median (nutrition-high)` = f3(`Median (nutrition-high)`),
                `Median (nutrition-low)` = f3(`Median (nutrition-low)`),
                Difference = f3(Difference),
                P = fmt_p(P), FDR = fmt_p(FDR))]
fwrite(t5out, file.path(TDIR, "Table4_immune_diff_by_nutrition_axis.csv"))

# ---------------- Table S3 按综合风险分分组（附循环论证说明） ----------------
r3 <- fread(file.path(PATHA$results, "composite_group_diff.tsv"))
setorder(r3, -diff)
ts3 <- r3[, .(immune_feature,
              `Median (high score)` = f3(median_high),
              `Median (low score)` = f3(median_low),
              Difference = f3(diff), P = fmt_p(p), FDR = fmt_p(fdr))]
fwrite(ts3, file.path(TDIR, "TableS3_immune_diff_by_composite_score.csv"))

# ---------------- Table S1 终点敏感性 ----------------
es <- fread(file.path(PATHA$results, "endpoint_sensitivity.tsv"))
ts1 <- es[, .(Endpoint = endpoint, Cohort = cohort, n, Events = events,
              `C-index` = f3(cindex), P = fmt_p(logrank_p),
              `Axes selected` = n_axis_sel, Note = note)]
fwrite(ts1, file.path(TDIR, "TableS1_endpoint_sensitivity.csv"))

# ---------------- Table S2 特征方案敏感性 ----------------
gs <- fread(file.path(PATHA$results, "geneset_sensitivity.tsv"))
ts2 <- gs[, .(Endpoint = endpoint, Training = train, Validation = test,
              Scheme = scheme, `n features` = n_feature, `n selected` = n_sel,
              `C train` = f3(c_train), `C test` = f3(c_test), P = fmt_p(p_test))]
fwrite(ts2, file.path(TDIR, "TableS2_geneset_sensitivity.csv"))

message("\n[09] 表格已写出 -> results/tables/")
print(list.files(TDIR))
message("\n[09] Table 1 预览："); print(t1)
message("\n[09] Table 4 预览："); print(t4)

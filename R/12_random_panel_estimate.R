# =====================================================================
# 12_random_panel_estimate.R
#   两项「给审稿人看的」补充分析，先出数再决定是否进正文：
#
#   A. 随机面板对照（random panel control）
#      - N1 表达匹配随机面板：对每条轴的每个基因，从「平均表达同一分位箱」中
#        随机抽 1 个基因替换，保持面板 16 个基因、且整体表达谱与真轴可比；
#      - N2 无约束随机面板：从全基因集均匀抽 16 个基因；
#      两种零分布各 K 个面板，用与主分析完全相同的 GSVA ssGSEA（同一次调用，
#      同参数）打分，再比较：
#        (i)  与 13 个 MSigDB 免疫特征的 Spearman 相关（均值 |rho| 与最大 |rho|）；
#        (ii) LUAD/LUSC 单变量 Cox（每 SD 的 logHR）与 Harrell C-index。
#      实证 P = (1 + #{null 不弱于 obs}) / (K + 1)。
#
#   B. ESTIMATE 纯度校正
#      - 用 ESTIMATE 作者站点（MD Anderson bioinformatics）的主表
#        rawdata/ESTIMATE_scores.txt 中 RNAseqV2 + 两种肺癌的行：Stromal /
#        Immune / ESTIMATE 三分量 + **tumor.purity**（作者给出的纯度），按 15 位
#        条码匹配；数值原样使用，不重算、不换算；
#      - 主校正协变量 = tumor.purity；另附 ESTIMATEScore（纯度单调反相关）与
#        StromalScore 单独校正，作为敏感性分析；
#      - 偏 Spearman：r_xy|z = (r_xy - r_xz r_yz) / sqrt((1-r_xz^2)(1-r_yz^2))，
#        P 用 t 分布 df = n - 3；raw 与 adjusted 必须在同一样本子集上计算，
#        否则 Δrho 无意义。
#
# 输出（均在 results/）：
#   random_panel_null_immune.tsv      轴 × 零模型 × 免疫特征：obs / null / p_emp
#   random_panel_null_survival.tsv    轴 × 零模型 × 队列：logHR、C-index 与零分布
#   random_panel_null_summary.tsv     汇总 4 项指标
#   random_panel_draws.tsv.gz         每次抽样的原始指标（可复算）
#   estimate_scores_matched.tsv       样本级 ESTIMATE 分数（含匹配标记）
#   estimate_axis_correlation.tsv     轴 × ESTIMATE 三分量
#   estimate_partial_correlation.tsv  轴 × 免疫特征：raw vs 校正
#   estimate_partial_summary.tsv      按「是否共享基因」分层的 Δrho 汇总
# 图：figures/_intermediate/figS_random_null.png、figS_purity_adjust.png
#
# 环境变量：PA_RANDOM_K 覆盖抽样数（默认 1000，冒烟测试用 20）
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

K <- as.integer(Sys.getenv("PA_RANDOM_K", "1000"))
SEED <- PARS$seed
set.seed(SEED)
message(sprintf("[12] 随机面板抽样数 K = %d；seed = %d", K, SEED))

d  <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
ex <- d$expr_log
axes <- names(d$genesets)

# ---------------- 0. 免疫特征（与 07 完全一致） ----------------
read_gmt <- function(path, want) {
  ln <- readLines(path)
  out <- list()
  for (l in ln) {
    p <- strsplit(l, "\t")[[1]]
    if (!(p[1] %in% want)) next
    out[[p[1]]] <- p[-c(1, 2)]
  }
  miss <- setdiff(want, names(out))
  if (length(miss)) stop("gmt 未找到: ", paste(miss, collapse = ", "))
  out
}
HALL <- file.path(PATHA$dir, "data/genesets/h.all.symbols.gmt")
GOBP <- file.path(PATHA$dir, "data/genesets/c5.go.bp.symbols.gmt")
IMM_GO <- c("GOBP_T_CELL_ACTIVATION", "GOBP_B_CELL_ACTIVATION",
            "GOBP_MACROPHAGE_ACTIVATION", "GOBP_NATURAL_KILLER_CELL_ACTIVATION",
            "GOBP_MYELOID_DENDRITIC_CELL_ACTIVATION",
            "GOBP_REGULATORY_T_CELL_DIFFERENTIATION",
            "GOBP_MAST_CELL_ACTIVATION", "GOBP_MONOCYTE_ACTIVATION",
            "GOBP_NEUTROPHIL_DEGRANULATION",
            "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION")
IMM_HALL <- c("HALLMARK_INFLAMMATORY_RESPONSE",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE",
              "HALLMARK_ALLOGRAFT_REJECTION")
imm_sets <- c(read_gmt(GOBP, IMM_GO), read_gmt(HALL, IMM_HALL))
imm_sets <- lapply(imm_sets, function(g) intersect(g, rownames(ex)))
imm_sets <- imm_sets[vapply(imm_sets, length, 0L) >= 10]

IMM_CACHE <- file.path(PATHA$work, "pa_scores_immune.rds")
if (file.exists(IMM_CACHE)) {
  imm <- readRDS(IMM_CACHE)
} else {
  imm <- ssgsea_score(ex, imm_sets)
  saveRDS(imm, IMM_CACHE)
}
imm <- as.data.frame(imm[rownames(imm) %in% colnames(ex), , drop = FALSE])
message(sprintf("[12] 免疫特征 %d 条；样本 %d", ncol(imm), nrow(imm)))

# ---------------- 1. 构造随机面板 ----------------
# 宇宙：方差非零的基因；轴基因本身从池中剔除（避免抽到自己）
v <- apply(ex, 1, stats::var)
pool_all <- names(v)[is.finite(v) & v > 0]
axis_genes <- unique(unlist(d$genesets))
pool <- setdiff(pool_all, axis_genes)
message(sprintf("[12] 抽样池 %d 基因（剔除轴基因 %d 个后）", length(pool), length(axis_genes)))

# 分位箱必须建在**全部基因**上：轴基因已被排除出 pool，若只在 pool 上建箱，
# 模板基因取不到箱号（NA），抽样会退化到全库 setdiff，复杂度从 O(2.5k) 升到 O(52k)
mu <- rowMeans(ex)
brk <- unique(stats::quantile(mu, probs = seq(0, 1, length.out = 21), na.rm = TRUE))
binc <- cut(mu, brk, include.lowest = TRUE, labels = FALSE)
names(binc) <- names(mu)
bin_list <- split(names(binc), binc)

make_matched <- function(template) {
  tb <- binc[match(template, names(binc))]
  out <- character(0)
  for (i in seq_along(template)) {
    cand <- setdiff(bin_list[[as.character(tb[i])]], c(axis_genes, out))
    if (!length(cand)) {
      # 极端情形：该表达箱内除模板基因外无可用基因 —— 退到相邻箱，
      # 仍保持「表达可比」，并在日志中可见（不静默改成全库随机）
      for (d in c(-1, 1, -2, 2, -3, 3)) {
        nb <- as.character(tb[i] + d)
        cand <- setdiff(bin_list[[nb]], c(axis_genes, out))
        if (length(cand)) break
      }
    }
    if (!length(cand)) cand <- setdiff(pool, c(axis_genes, out))
    out <- c(out, sample(cand, 1))
  }
  out
}

# 预分配，避免在 6,000 次循环里逐步 rbind / 追加（O(n^2) 复制）
n_pan <- length(axes) * K * 2L
panels <- vector("list", n_pan)
pid <- character(n_pan); pax <- character(n_pan)
pty <- character(n_pan); pdr <- integer(n_pan)
ii <- 0L
for (a in axes) {
  tg <- d$genesets[[a]]
  for (k in seq_len(K)) {
    ii <- ii + 1L
    pid[ii] <- sprintf("N1_%s_%04d", a, k); pax[ii] <- a
    pty[ii] <- "expr_matched"; pdr[ii] <- k
    panels[[ii]] <- make_matched(tg)
    ii <- ii + 1L
    pid[ii] <- sprintf("N2_%s_%04d", a, k); pax[ii] <- a
    pty[ii] <- "unconstrained"; pdr[ii] <- k
    panels[[ii]] <- sample(pool, length(tg))
  }
}
names(panels) <- pid
meta_p <- data.table(panel_id = pid, axis = pax, null_type = pty, draw = pdr)
# 真轴也一起打分，保证「同一次 GSVA 调用、同一参数」，并作为 QC 复现主分析分数
obs_sets <- lapply(axes, function(a) d$genesets[[a]])
names(obs_sets) <- paste0("OBS_", axes)
all_sets <- c(obs_sets, panels)
message(sprintf("[12] 待打分基因集 %d 个（含 3 条真轴）", length(all_sets)))

# 中间结果缓存：6,003 个基因集的 ssGSEA + 两个队列的 Cox/C-index 合计约 18 分钟，
# 重复调试时不重算（缓存键含 K 与 seed，改参数自动失效）
CACHE <- file.path(PATHA$work, sprintf("pa_random_panels_K%d_seed%d.rds", K, SEED))
if (file.exists(CACHE)) {
  message("[12] 读取中间缓存：", basename(CACHE))
  .c <- readRDS(CACHE)
  ps <- .c$ps; mL <- .c$mL; mV <- .c$mV
} else {
  t0 <- Sys.time()
  ps <- ssgsea_score(ex, all_sets)
  message(sprintf("[12] ssGSEA 完成，用时 %.1f 分钟", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# QC：同参数下真轴分数应与主分析 pa_scores.rds 完全一致
qc <- vapply(axes, function(a) {
  stats::cor(ps[[paste0("OBS_", a)]], sc[[a]], method = "spearman")
}, 0)
message("[12] QC 真轴复算 Spearman（应≈1）: ",
        paste(sprintf("%s=%.4f", names(qc), qc), collapse = ", "))
if (any(qc < 0.999)) warning("真轴复算与主分析不一致，检查 ssGSEA 参数")

# ---------------- 2. 指标计算 ----------------
# 2a 免疫特征相关：Spearman = 秩上的 Pearson，一次矩阵乘法出全部面板
cm <- intersect(rownames(imm), colnames(ex))
I <- as.matrix(imm[cm, , drop = FALSE])
P <- as.matrix(ps[cm, , drop = FALSE])
Ri <- scale(apply(I, 2, rank))
Rp <- scale(apply(P, 2, rank))
rho <- crossprod(Rp, Ri) / (nrow(I) - 1)   # 面板 × 特征
message(sprintf("[12] 相关矩阵 %d 面板 × %d 特征（n=%d）", nrow(rho), ncol(rho), nrow(I)))

# 2b 生存指标（LUAD 训练 / LUSC 验证）
surv_metric <- function(dt, mat) {
  s <- intersect(dt$sample, rownames(mat))
  dt <- dt[match(s, dt$sample)]
  m <- mat[s, , drop = FALSE]
  m <- scale(m)   # 每 SD
  out <- vapply(seq_len(ncol(m)), function(j) {
    x <- m[, j]
    f <- suppressWarnings(survival::coxph(survival::Surv(dt$time, dt$status) ~ x))
    b <- summary(f)$coefficients[1, c(1, 5)]
    c(logHR = b[1], p = b[2],
      C = cindex(dt$time, dt$status, as.numeric(x)))
  }, c(logHR = 0, p = 0, C = 0))
  list(n = length(s), n_ev = sum(dt$status == 1), m = t(out))
}
tr <- merge(d$surv_luad, sc, by = "sample")
va <- merge(d$surv_lusc, sc, by = "sample")
message(sprintf("[12] 队列：LUAD %d 例（事件 %d）；LUSC %d 例（事件 %d）",
                nrow(tr), sum(tr$status == 1), nrow(va), sum(va$status == 1)))
if (!exists("mL")) {
  t1 <- Sys.time()
  mL <- surv_metric(tr, P)
  mV <- surv_metric(va, P)
  message(sprintf("[12] 生存指标完成（%.1f 分钟）：LUAD n=%d（事件 %d）；LUSC n=%d（事件 %d）",
                  as.numeric(difftime(Sys.time(), t1, units = "mins")),
                  mL$n, mL$n_ev, mV$n, mV$n_ev))
  saveRDS(list(ps = ps, mL = mL, mV = mV), CACHE)
  message("[12] 中间结果已缓存：", basename(CACHE))
}

# 表达匹配 QC：面板基因的平均表达分布 vs 模板轴的分布（排序后 Spearman；
# 观察到的轴与自身比即为 1，用于确认匹配确实生效）
set_axis <- c(setNames(axes, paste0("OBS_", axes)),
              setNames(meta_p$axis, meta_p$panel_id))
muex <- rowMeans(ex)      # 只算一次；在 6,003 次迭代里重算是 O(6000 x 60k x 1011)
expr_rho <- vapply(names(all_sets), function(id) {
  tg <- all_sets[[id]]
  ref <- d$genesets[[set_axis[[id]]]]
  stats::cor(sort(muex[tg]), sort(muex[ref]), method = "spearman")
}, 0)

# ---------------- 3. 零分布比较 ----------------
draws <- data.table(panel_id = names(all_sets), axis = NA_character_,
                    null_type = "observed", draw = 0L)
draws[grepl("^N[12]_", panel_id), `:=`(
  axis = sub("^N[12]_([a-z_]+)_[0-9]+$", "\\1", panel_id),
  null_type = ifelse(grepl("^N1_", panel_id), "expr_matched", "unconstrained"),
  draw = as.integer(sub("^.*_([0-9]+)$", "\\1", panel_id)))]
draws[grepl("^OBS_", panel_id), axis := sub("^OBS_", "", panel_id)]
draws[, `:=`(rho_mean = rowMeans(rho),
             rho_absmean = rowMeans(abs(rho)),
             rho_absmax  = apply(abs(rho), 1, max),
             expr_match_rho = expr_rho[panel_id],
             logHR_luad = mL$m[, "logHR"], C_luad = mL$m[, "C"],
             p_luad = mL$m[, "p"],
             logHR_lusc = mV$m[, "logHR"], C_lusc = mV$m[, "C"],
             p_lusc = mV$m[, "p"])]
fwrite(draws, file.path(PATHA$results, "random_panel_draws.tsv.gz"), sep = "\t")

emp_p <- function(null, obs, two_sided = FALSE) {
  if (two_sided) (1 + sum(abs(null) >= abs(obs))) / (length(null) + 1)
  else (1 + sum(null >= obs)) / (length(null) + 1)
}

summ <- rbindlist(lapply(axes, function(a) {
  rbindlist(lapply(c("expr_matched", "unconstrained"), function(nt) {
    nd <- draws[axis == a & null_type == nt]
    ob <- draws[panel_id == paste0("OBS_", a)]
    rbindlist(list(
      data.table(metric = "immune_rho_mean", obs = ob$rho_mean,
                 null_mean = mean(nd$rho_mean), null_sd = sd(nd$rho_mean),
                 q025 = quantile(nd$rho_mean, .025), q975 = quantile(nd$rho_mean, .975),
                 z = (ob$rho_mean - mean(nd$rho_mean)) / sd(nd$rho_mean),
                 p_emp = emp_p(nd$rho_mean, ob$rho_mean), n_null = nrow(nd)),
      data.table(metric = "immune_rho_absmean", obs = ob$rho_absmean,
                 null_mean = mean(nd$rho_absmean), null_sd = sd(nd$rho_absmean),
                 q025 = quantile(nd$rho_absmean, .025), q975 = quantile(nd$rho_absmean, .975),
                 z = (ob$rho_absmean - mean(nd$rho_absmean)) / sd(nd$rho_absmean),
                 p_emp = emp_p(nd$rho_absmean, ob$rho_absmean), n_null = nrow(nd)),
      data.table(metric = "immune_rho_absmax", obs = ob$rho_absmax,
                 null_mean = mean(nd$rho_absmax), null_sd = sd(nd$rho_absmax),
                 q025 = quantile(nd$rho_absmax, .025), q975 = quantile(nd$rho_absmax, .975),
                 z = (ob$rho_absmax - mean(nd$rho_absmax)) / sd(nd$rho_absmax),
                 p_emp = emp_p(nd$rho_absmax, ob$rho_absmax), n_null = nrow(nd)),
      data.table(metric = "os_Cindex_luad", obs = ob$C_luad,
                 null_mean = mean(nd$C_luad), null_sd = sd(nd$C_luad),
                 q025 = quantile(nd$C_luad, .025), q975 = quantile(nd$C_luad, .975),
                 z = (ob$C_luad - mean(nd$C_luad)) / sd(nd$C_luad),
                 p_emp = emp_p(nd$C_luad, ob$C_luad, two_sided = TRUE), n_null = nrow(nd))
    ))[, `:=`(axis = a, null_type = nt)]
  }))
}))
setcolorder(summ, c("axis", "null_type", "metric"))
fwrite(summ, file.path(PATHA$results, "random_panel_null_summary.tsv"), sep = "\t")

# 逐免疫特征的 obs vs 零分布（附基因重叠度，沿用 07 的分层口径）
imm_null <- rbindlist(lapply(axes, function(a) {
  rbindlist(lapply(c("expr_matched", "unconstrained"), function(nt) {
    ids_n <- draws[axis == a & null_type == nt, panel_id]
    ids_o <- paste0("OBS_", a)
    rbindlist(lapply(colnames(rho), function(f) {
      ob <- rho[ids_o, f]
      nd <- rho[ids_n, f]
      data.table(axis = a, null_type = nt, immune_feature = f,
                 overlap_genes = length(intersect(d$genesets[[a]], imm_sets[[f]])),
                 obs_rho = ob, null_mean = mean(nd), null_sd = sd(nd),
                 z = (ob - mean(nd)) / sd(nd),
                 p_emp = emp_p(nd, ob))
    }))
  }))
}))
imm_null[, fdr := p.adjust(p_emp, "BH")]
fwrite(imm_null, file.path(PATHA$results, "random_panel_null_immune.tsv"), sep = "\t")

# 生存：obs 与零分布并列（两队列）
surv_null <- rbindlist(lapply(axes, function(a) {
  rbindlist(lapply(c("expr_matched", "unconstrained"), function(nt) {
    nd <- draws[axis == a & null_type == nt]
    ob <- draws[panel_id == paste0("OBS_", a)]
    rbindlist(list(
      data.table(cohort = "LUAD", metric = "logHR", obs = ob$logHR_luad,
                 null_mean = mean(nd$logHR_luad), null_sd = sd(nd$logHR_luad),
                 q025 = quantile(nd$logHR_luad, .025), q975 = quantile(nd$logHR_luad, .975),
                 p_emp = emp_p(nd$logHR_luad, ob$logHR_luad, two_sided = TRUE)),
      data.table(cohort = "LUSC", metric = "logHR", obs = ob$logHR_lusc,
                 null_mean = mean(nd$logHR_lusc), null_sd = sd(nd$logHR_lusc),
                 q025 = quantile(nd$logHR_lusc, .025), q975 = quantile(nd$logHR_lusc, .975),
                 p_emp = emp_p(nd$logHR_lusc, ob$logHR_lusc, two_sided = TRUE))
    ))[, `:=`(axis = a, null_type = nt)]
  }))
}))
surv_null[, fdr := p.adjust(p_emp, "BH")]
fwrite(surv_null, file.path(PATHA$results, "random_panel_null_survival.tsv"), sep = "\t")

message("\n[12] A. 随机面板对照汇总：")
print(summ[, .(axis, null_type, metric, obs = round(obs, 3),
               null_mean = round(null_mean, 3),
               null_95 = sprintf("%.3f-%.3f", q025, q975),
               z = round(z, 2), p_emp = signif(p_emp, 3))])

# ---------------- 4. ESTIMATE ----------------
# 用 ESTIMATE 作者站点（bioinformatics.mdanderson.org/estimate）发布的
# rawdata/ESTIMATE_scores.txt 主表：每行含 Platform / Disease / 三分量 / tumor.purity。
# 只取 RNAseqV2 + 两种肺癌，按 15 位条码匹配；不改写任何数值。
EST_FILE <- file.path(PATHA$dir, "data/external/estimate/ESTIMATE_scores.txt")
EST <- fread(EST_FILE)
est <- EST[Platform == "RNAseqV2" &
             Disease %in% c("lung adenocarcinoma", "lung squamous cell carcinoma"),
           .(barcode15 = substr(ID, 1, 15),
             Stromal = Stromal_score, Immune = Immune_score,
             ESTIMATE = ESTIMATE_score, purity = tumor.purity,
             src_project = fifelse(Disease == "lung adenocarcinoma", "TCGA-LUAD", "TCGA-LUSC"))]
setorder(est, barcode15)

sm <- data.table(sample = rownames(ps))
sm[, barcode15 := substr(sample, 1, 15)]
# 队列归属以生存表为准（与主分析 500/491 完全一致）
sm[, cohort := fifelse(sample %in% tr$sample, "LUAD",
                fifelse(sample %in% va$sample, "LUSC", "no_OS"))]
dup <- est$barcode15[duplicated(est$barcode15)]
if (length(dup)) stop("ESTIMATE 表存在重复条码: ", paste(head(dup, 3), collapse = ", "))
sm <- merge(sm, est[, .(barcode15, Stromal, Immune, ESTIMATE, purity)], by = "barcode15", all.x = TRUE)
sm[, matched := !is.na(ESTIMATE)]
fwrite(sm, file.path(PATHA$results, "estimate_scores_matched.tsv"), sep = "\t")
message(sprintf("\n[12] B. ESTIMATE 匹配：LUAD %d/%d，LUSC %d/%d",
                sm[cohort == "LUAD" & matched == TRUE, .N], nrow(tr),
                sm[cohort == "LUSC" & matched == TRUE, .N], nrow(va)))

# 4a 轴 × ESTIMATE 三分量（轴分数取同一次 ssGSEA 的 OBS_ 列，保证与零分布可比）
obs_cols <- paste0("OBS_", axes)
S <- as.data.frame(ps[sm$sample, obs_cols, drop = FALSE], check.names = FALSE)
names(S) <- axes; rownames(S) <- sm$sample
ok <- sm$matched
r_est <- rbindlist(lapply(axes, function(a) {
  rbindlist(lapply(c("Stromal", "Immune", "ESTIMATE"), function(k) {
    ct <- suppressWarnings(stats::cor.test(S[ok, a], sm[ok][[k]], method = "spearman"))
    data.table(axis = a, score_type = k, rho = unname(ct$estimate),
               p = ct$p.value, n = sum(ok))
  }))
}))
fwrite(r_est, file.path(PATHA$results, "estimate_axis_correlation.tsv"), sep = "\t")

# 4b 偏 Spearman：轴 × 免疫特征，控制 ESTIMATEScore（敏感性：仅 StromalScore）
psp <- function(x, y, z) {
  rxy <- stats::cor(x, y, method = "spearman")
  rxz <- stats::cor(x, z, method = "spearman")
  ryz <- stats::cor(y, z, method = "spearman")
  num <- rxy - rxz * ryz
  den <- sqrt((1 - rxz^2) * (1 - ryz^2))
  r <- num / den
  r[abs(den) < 1e-12] <- NA_real_
  df <- length(x) - 3
  t <- r * sqrt(df / (1 - r^2))
  list(r = r, p = 2 * stats::pt(-abs(t), df))
}
sub <- sm[matched == TRUE]
S2 <- as.data.frame(ps[sub$sample, obs_cols, drop = FALSE], check.names = FALSE)
names(S2) <- axes; rownames(S2) <- sub$sample
I2 <- as.data.frame(imm[sub$sample, , drop = FALSE])

part <- rbindlist(lapply(axes, function(a) {
  rbindlist(lapply(colnames(I2), function(f) {
    raw <- suppressWarnings(stats::cor.test(S2[[a]], I2[[f]], method = "spearman"))
    pe <- psp(S2[[a]], I2[[f]], sub$ESTIMATE)
    pp <- psp(S2[[a]], I2[[f]], sub$purity)
    ps2 <- psp(S2[[a]], I2[[f]], sub$Stromal)
    data.table(axis = a, immune_feature = f,
               overlap_genes = length(intersect(d$genesets[[a]], imm_sets[[f]])),
               n = nrow(sub),
               rho_raw = unname(raw$estimate), p_raw = raw$p.value,
               rho_adj_estimate = pe$r, p_adj_estimate = pe$p,
               rho_adj_purity = pp$r, p_adj_purity = pp$p,
               rho_adj_stromal = ps2$r, p_adj_stromal = ps2$p)
  }))
}))
part[, `:=`(fdr_raw = p.adjust(p_raw, "BH"),
            fdr_adj_estimate = p.adjust(p_adj_estimate, "BH"),
            fdr_adj_purity = p.adjust(p_adj_purity, "BH"),
            fdr_adj_stromal = p.adjust(p_adj_stromal, "BH"),
            delta_purity = rho_raw - rho_adj_purity,
            delta_estimate = rho_raw - rho_adj_estimate,
            delta_stromal = rho_raw - rho_adj_stromal,
            kind = ifelse(overlap_genes > 0, "homologous (shared genes)",
                          "independent (no shared genes)"))]
fwrite(part, file.path(PATHA$results, "estimate_partial_correlation.tsv"), sep = "\t")

psum <- part[, .(n_pairs = .N,
                 rho_raw_median = median(rho_raw),
                 rho_adj_purity_median = median(rho_adj_purity),
                 delta_purity_median = median(delta_purity),
                 rho_adj_est_median = median(rho_adj_estimate),
                 delta_est_median = median(delta_estimate),
                 rho_adj_str_median = median(rho_adj_stromal),
                 delta_str_median = median(delta_stromal),
                 n_sig_raw = sum(fdr_raw < 0.05),
                 n_sig_adj_purity = sum(fdr_adj_purity < 0.05)),
             by = .(axis, kind)]
fwrite(psum, file.path(PATHA$results, "estimate_partial_summary.tsv"), sep = "\t")

message("\n[12] 轴 × ESTIMATE（Spearman）：")
print(r_est[, .(axis, score_type, rho = round(rho, 3), p = signif(p, 2))])
message("\n[12] raw vs 纯度校正（按是否共享基因分层）：")
print(psum[, .(axis, kind, n_pairs,
               raw = round(rho_raw_median, 3),
               adj_purity = round(rho_adj_purity_median, 3),
               delta = round(delta_purity_median, 3),
               adj_est = round(rho_adj_est_median, 3),
               adj_str = round(rho_adj_str_median, 3),
               sig_raw = n_sig_raw, sig_adj = n_sig_adj_purity)])

# ---------------- 5. 图（供人工判断，暂不进投稿包） ----------------
suppressPackageStartupMessages({ library(ggplot2) })
INT <- file.path(PATHA$figures, "_intermediate")
dir.create(INT, showWarnings = FALSE, recursive = TRUE)

pd <- rbind(
  draws[null_type != "observed", .(axis, null_type, value = rho_absmean,
                                   set = "null")],
  draws[null_type == "observed", .(axis, null_type = "observed",
                                   value = rho_absmean, set = "observed")]
)
p1 <- ggplot(pd, aes(x = value, fill = null_type)) +
  geom_histogram(bins = 40, alpha = .6, position = "identity") +
  geom_vline(data = pd[set == "observed"], aes(xintercept = value),
             colour = "#C0392B", linewidth = .9) +
  facet_wrap(~axis, scales = "free_y") +
  labs(x = "mean |Spearman rho| across 13 immune features", y = "panels",
       title = "Random-panel null vs observed axes") +
  theme_bw(base_size = 11)
ggsave(file.path(INT, "figS_random_null.png"), p1, width = 8, height = 4.2, dpi = 150)

p2 <- ggplot(part, aes(x = rho_raw, y = rho_adj_purity, colour = kind)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 1.6) +
  facet_wrap(~axis) +
  labs(x = "Spearman rho (raw)", y = "partial rho (adjusted for ESTIMATE tumour purity)",
       title = "Purity adjustment: raw vs ESTIMATE-purity-adjusted") +
  theme_bw(base_size = 11)
ggsave(file.path(INT, "figS_purity_adjust.png"), p2, width = 8, height = 3.4, dpi = 150)

message("\n[12] 完成 -> results/random_panel_*.tsv、estimate_*.tsv；figures/_intermediate/figS_*.png")

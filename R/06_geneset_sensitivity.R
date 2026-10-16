# =====================================================================
# 06_geneset_sensitivity.R — 两个变量的实测对比：
#   (A) 特征方案：S0 先验 3 轴(16 基因) / S1 MSigDB 6 条真实通路 / S2 48 基因单基因级
#   (B) 训练方向：LUAD→LUSC 与 LUSC→LUAD 双向（质疑「LUAD 当训练集」这一默认）
# 终点：OS 与 PFI 各跑一遍。
# 输出：results/geneset_sensitivity.tsv, results/axis_direction_double.tsv
# 说明：S1 使用的 6 条通路均为真实 MSigDB 条目（名称已逐条 grep 核实存在），
#       不虚构通路名；S2 只用原始 48 个先验基因，不做任何额外挑选。
#       【版本口径·2026-10-16 已确证】本地 gmt 为 MSigDB **v2024.1.Hs**。
#       确证方式：与 Broad **公开**归档（无需登录）
#         https://data.broadinstitute.org/gsea-msigdb/msigdb/release/<版本>/<集合>.v<版本>.symbols.gmt
#       做 md5 逐字节比对——本地 c5.go.bp（fd535c44…，7,608 集）与 h.all（d3606c32…，50 集）
#       均**唯一**命中 2024.1.Hs；六个候选版本 v7.5.1 / 2023.1.Hs / 2023.2.Hs /
#       2024.1.Hs / 2025.1.Hs / 2026.1.Hs 中不存在第二个命中。
#       复跑：bash tools/verify_msigdb_version.sh → data/_vercheck/msigdb_version_verdict.txt
#       （更正：早期「无法由文件内容确证、归档需登录」的判断为**误判**，已作废；详见
#        manuscript/02_methods.md 附录 C 第 2.5 行。）
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

d   <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc3 <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
ex  <- d$expr_log
axes <- names(d$genesets)

# ---- gmt 读取 ----
read_gmt <- function(path, want) {
  ln <- readLines(path)
  out <- list()
  for (l in ln) {
    p <- strsplit(l, "\t")[[1]]
    if (!(p[1] %in% want)) next
    out[[p[1]]] <- p[-c(1, 2)]
  }
  missing <- setdiff(want, names(out))
  if (length(missing)) stop("gmt 中未找到: ", paste(missing, collapse = ", "))
  out
}
HALL <- file.path(PATHA$dir, "data/genesets/h.all.symbols.gmt")
GOBP <- file.path(PATHA$dir, "data/genesets/c5.go.bp.symbols.gmt")

gs_msig <- c(
  read_gmt(HALL, c("HALLMARK_COAGULATION", "HALLMARK_INFLAMMATORY_RESPONSE",
                   "HALLMARK_ALLOGRAFT_REJECTION")),
  read_gmt(GOBP, c("GOBP_ACUTE_PHASE_RESPONSE", "GOBP_NEUTROPHIL_DEGRANULATION",
                   "GOBP_ADAPTIVE_IMMUNE_RESPONSE"))
)
gs_msig <- lapply(gs_msig, function(g) intersect(g, rownames(ex)))
message("[06] S1 通路规模: ",
        paste(sprintf("%s=%d", sub("^(HALLMARK|GOBP)_", "", names(gs_msig)),
                      vapply(gs_msig, length, 0L)), collapse = ", "))

# ---- 特征矩阵（行 = 样本）----
S0 <- as.data.frame(sc3[, axes]); rownames(S0) <- sc3$sample
S1_CACHE <- file.path(PATHA$work, "pa_scores_msig6.rds")
if (file.exists(S1_CACHE)) {
  S1 <- readRDS(S1_CACHE)
  message("[06] S1 打分读缓存")
} else {
  S1 <- ssgsea_score(ex, gs_msig)
  saveRDS(S1, S1_CACHE)
}
g48 <- intersect(unique(unlist(d$genesets)), rownames(ex))
S2 <- as.data.frame(t(ex[g48, , drop = FALSE]))
FEATS <- list(S0_prior16 = S0, S1_msigdb6 = S1, S2_genes48 = S2)

# ---- 评估 ----
eval_scheme <- function(feat, sv_tr, sv_va, endpoint, scheme) {
  getcol <- function(sv) {
    sv[, time := get(paste0(endpoint, "_time"))]
    sv[, status := get(endpoint)]
    sv[!is.na(time) & !is.na(status) & time > 0]
  }
  tr <- getcol(copy(sv_tr)); va <- getcol(copy(sv_va))
  tr <- tr[sample %in% rownames(feat)]; va <- va[sample %in% rownames(feat)]
  ftr <- as.matrix(feat[tr$sample, , drop = FALSE])
  fva <- as.matrix(feat[va$sample, , drop = FALSE])
  mu <- colMeans(ftr); sdv <- apply(ftr, 2, stats::sd)
  keep <- is.finite(sdv) & sdv > 0
  ftr <- scale(ftr[, keep, drop = FALSE], center = mu[keep], scale = sdv[keep])
  fva <- scale(fva[, keep, drop = FALSE], center = mu[keep], scale = sdv[keep])
  fit <- suppressWarnings(fit_lasso_cox(as.data.frame(ftr), tr$time, tr$status,
                                        nfolds = PARS$lasso_nfolds, seed = PARS$seed))
  rs_tr <- risk_score(as.data.frame(ftr), fit$coef)
  rs_va <- risk_score(as.data.frame(fva), fit$coef)
  nz <- sum(fit$coef != 0)
  ctr <- if (nz && stats::sd(rs_tr) > 0) cindex(tr$time, tr$status, rs_tr) else 0.5
  cva <- if (nz && stats::sd(rs_va) > 0) cindex(va$time, va$status, rs_va) else 0.5
  pva <- NA_real_
  if (nz && stats::sd(rs_va) > 0) {
    grp <- ifelse(rs_va > median(rs_tr), "High", "Low")
    if (length(unique(grp)) < 2) grp <- ifelse(rs_va > median(rs_va), "High", "Low")
    if (length(unique(grp)) >= 2) pva <- logrank(va$time, va$status, grp)
  }
  data.table(scheme = scheme, endpoint = endpoint,
             n_feature = ncol(ftr), n_sel = nz,
             n_train = nrow(tr), n_test = nrow(va),
             c_train = round(ctr, 4), c_test = round(cva, 4),
             p_test = signif(pva, 3),
             s_used = ifelse(is.character(fit$s_used), fit$s_used,
                             sprintf("%.4g", fit$s_used)))
}

dirs <- list(c("LUAD", "LUSC"), c("LUSC", "LUAD"))
rows <- list()
for (ep in c("OS", "PFI")) {
  for (dd in dirs) {
    sv_tr <- if (dd[1] == "LUAD") d$surv_luad else d$surv_lusc
    sv_va <- if (dd[2] == "LUAD") d$surv_luad else d$surv_lusc
    for (nm in names(FEATS)) {
      r <- eval_scheme(FEATS[[nm]], sv_tr, sv_va, ep, nm)
      r[, `:=`(train = dd[1], test = dd[2])]
      rows[[paste(ep, dd[1], nm)]] <- r
      message(sprintf("  %s %s->%s %-12s 特征%3d 选中%2d  C_train=%.3f  C_test=%.3f  P_test=%.3g",
                      ep, dd[1], dd[2], nm, r$n_feature, r$n_sel,
                      r$c_train, r$c_test, r$p_test))
    }
  }
}
out <- rbindlist(rows)
setcolorder(out, c("endpoint", "train", "test", "scheme", "n_feature", "n_sel",
                   "n_train", "n_test", "c_train", "c_test", "p_test", "s_used"))
fwrite(out, file.path(PATHA$results, "geneset_sensitivity.tsv"), sep = "\t")
message("\n[06] 完成 -> results/geneset_sensitivity.tsv")
print(out)

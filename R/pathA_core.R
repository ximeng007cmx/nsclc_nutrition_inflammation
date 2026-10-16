# =====================================================================
# pathA_core.R — 路径 A 核心函数库（数据读取 / ssGSEA / LASSO-Cox / 评估）
# 设计原则：
#   * 每个函数纯函数式、可单独自测；
#   * GSVA 2.x 用 param 对象式 API（ssgseaParam + gsva）；
#   * 不依赖未安装包（不用 survminer/estimate/WGCNA/survAUC）。
# =====================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------- 1. 数据读取 ----------

#' 读取 Toil probeMap（实为纯文本 TSV：id gene chrom chromStart chromEnd strand）
read_probemap <- function(path) {
  m <- fread(path, header = TRUE, sep = "\t")
  setnames(m, c("id", "gene"), c("ensg", "symbol"))
  m[, .(ensg, symbol)]
}

#' 读取子集表达矩阵（基因 × 样本，首列 gene = ENSG-versioned）
read_expr <- function(path) {
  dt <- fread(path, header = TRUE, sep = "\t")
  mat <- as.matrix(dt[, -1])
  rownames(mat) <- dt$gene
  mat
}

#' ENSG -> symbol，重复 symbol 取行均值
map_to_symbol <- function(mat, probemap) {
  hit <- probemap[symbol != "", .(ensg, symbol)]
  ix <- match(rownames(mat), hit$ensg)
  sym <- hit$symbol[ix]
  keep <- !is.na(sym)
  m2 <- mat[keep, , drop = FALSE]
  sym <- sym[keep]
  dt <- as.data.table(m2); dt[, symbol := sym]
  out <- dt[, lapply(.SD, mean), by = symbol]
  m3 <- as.matrix(out[, -1])
  rownames(m3) <- out$symbol
  m3
}

#' FPKM -> log2(FPKM + 1)
log2_fpkm <- function(m) log2(m + 1)

#' Toil 泛癌矩阵存储 log2(FPKM + 0.001)（下界 -9.9658 = log2(0.001)），
#' 先还原 FPKM 再取 log2(FPKM+1)，避免双重对数
toil_fix <- function(m) {
  if (min(m, na.rm = TRUE) < -5) {
    message("     检测到 Toil log2(FPKM+0.001) 编码，还原为 log2(FPKM+1)")
    fpkm <- 2^m - 0.001
    fpkm[fpkm < 0] <- 0
    m <- log2(fpkm + 1)
  }
  m
}

# ---------- 2. ssGSEA 打分 ----------

#' GSVA 2.x ssGSEA。返回 基因集 × 样本 分数矩阵。
ssgsea_score <- function(expr, genesets) {
  stopifnot(requireNamespace("GSVA", quietly = TRUE))
  p <- GSVA::ssgseaParam(exprData = as.matrix(expr), geneSets = genesets)
  s <- GSVA::gsva(p, verbose = FALSE)
  as.data.frame(t(s))  # 样本 × 基因集
}

# ---------- 3. 生存数据构造 ----------

#' 把临床表 + 样本元数据合成建模表（样本级）
#' endpoint: "OS" / "DSS" / "PFI"
build_surv_table <- function(meta, clinical, endpoint = "OS") {
  stopifnot(endpoint %in% c("OS", "DSS", "PFI"))
  meta <- as.data.table(meta)
  setnames(meta, c("sample", "project", "patient"))
  cl <- as.data.table(clinical)
  keep <- c("cancer_type", "patient", "gender", "stage", "vital_status",
            "OS_time", "OS", "DSS_time", "DSS", "PFI_time", "PFI", "age")
  keep <- intersect(keep, names(cl))
  cl <- cl[, ..keep]
  m <- merge(meta, cl, by.x = "patient", by.y = "patient", all.x = TRUE, sort = FALSE)
  m[, time := m[[paste0(endpoint, "_time")]]]
  m[, status := m[[endpoint]]]
  m <- m[!is.na(time) & !is.na(status)]
  m[, time := as.numeric(time)]
  m[, status := as.numeric(status)]
  m[time > 0]
}

# ---------- 4. LASSO-Cox ----------

#' 在训练集上：univariate Cox 初筛(可选) -> cv.glmnet(family="cox")
fit_lasso_cox <- function(scores, time, status, nfolds = 10, alpha = 1,
                          univariate_pmax = NULL, seed = 1) {
  set.seed(seed)
  x <- as.matrix(scores)
  y <- survival::Surv(time, status)

  keep <- seq_len(ncol(x))
  if (!is.null(univariate_pmax)) {
    pv <- vapply(keep, function(j) {
      summary(survival::coxph(y ~ x[, j]))$coefficients[1, "Pr(>|z|)"]
    }, 0)
    keep <- keep[pv <= univariate_pmax]
    if (length(keep) == 0) stop("univariate Cox 初筛后无变量入选")
    x <- x[, keep, drop = FALSE]
  }

  cvfit <- glmnet::cv.glmnet(x, y, family = "cox", alpha = alpha,
                             nfolds = nfolds, standardize = TRUE,
                             cox.ties = "efron")
  cf <- as.vector(stats::coef(cvfit, s = "lambda.min"))
  s_used <- "lambda.min"
  if (all(cf == 0)) {
    # 诚实回退：lambda.min 全零（信号弱）→ 取「最大的仍保留非零系数」的 lambda，
    # 并显式打标 fallback，供正文/结果如实报告
    bet <- cvfit$glmnet.fit$beta
    nzcol <- apply(as.matrix(bet), 2, function(cc) any(cc != 0))
    cand <- which(nzcol & cvfit$lambda <= cvfit$lambda.min)
    if (!length(cand)) cand <- which(nzcol)
    if (!length(cand)) stop("整个 lambda 路径均无非零系数，模型不可拟合")
    j <- max(cand)
    cf <- as.vector(bet[, j])
    s_used <- cvfit$lambda[j]
    warning(sprintf("lambda.min 系数全零（信号弱），回退 lambda=%.4g", s_used))
  }
  names(cf) <- colnames(x)
  list(cvfit = cvfit, coef = cf, kept = keep,
       lambda_min = cvfit$lambda.min, s_used = s_used)
}

#' 用既定系数给任意样本算风险分
risk_score <- function(new_scores, coef) {
  x <- as.matrix(new_scores)[, names(coef), drop = FALSE]
  as.vector(x %*% coef)
}

# ---------- 5. 评估 ----------

#' Harrell C-index（风险分语义：分数越高风险越高）
#' 注意：survival::concordance(Surv~x) 把 x 当保护变量编码（x 越大生存越长），
#' 实测风险信号返回 1-C（与 coxph summary 的 C 互补），故此处取 1-。
cindex <- function(time, status, score) {
  1 - survival::concordance(survival::Surv(time, status) ~ score)$concordance
}

#' C-index 及 bootstrap 95% CI
cindex_boot <- function(time, status, score, n = 500, seed = 1) {
  set.seed(seed)
  nobs <- length(time)
  est <- cindex(time, status, score)
  bs <- vapply(seq_len(n), function(i) {
    k <- sample.int(nobs, nobs, replace = TRUE)
    cindex(time[k], status[k], score[k])
  }, 0)
  ci <- stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE)
  c(est = est, lo = unname(ci[1]), hi = unname(ci[2]))
}

#' log-rank 检验（高/低风险组）
logrank <- function(time, status, group) {
  sd <- survival::survdiff(survival::Surv(time, status) ~ group)
  pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)
}

#' 多变量 Cox：风险分 + 临床协变量；返回整理好的 data.frame
cox_multivariate <- function(dt, covars = c("age", "gender", "stage_late")) {
  dt <- as.data.table(dt)
  rhs <- paste(c("risk_score", covars), collapse = " + ")
  f <- stats::as.formula(paste("survival::Surv(time, status) ~", rhs))
  fit <- survival::coxph(f, data = dt)
  s <- summary(fit)
  out <- data.table(
    term = rownames(s$coefficients),
    HR = s$conf.int[, "exp(coef)"],
    lo = s$conf.int[, "lower .95"],
    hi = s$conf.int[, "upper .95"],
    p = s$coefficients[, "Pr(>|z|)"]
  )
  list(fit = fit, table = out)
}

#' 临床分期二分（I/II 早期 vs III/IV 晚期）
stage_binary <- function(stage) {
  s <- toupper(gsub("[^A-Za-z0-9 ]", "", as.character(stage)))
  out <- rep(NA_character_, length(s))
  out[grepl("STAGE (I|II)(A|B|C)?$", s) | s %in% c("I", "II", "IA", "IB", "IIA", "IIB", "IIC")] <- "early"
  out[grepl("STAGE (III|IV)", s) | s %in% c("III", "IIIA", "IIIB", "IIIC", "IV")] <- "late"
  out
}

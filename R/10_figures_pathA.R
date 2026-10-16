# =====================================================================
# R/10_figures_pathA.R — 生成路线 1 的稿件图（Figure 2–5）
#
# 与 A1 的图片重编号方案对应：
#   Figure 1  研究流程       → 由 tools/make_fig1_workflow.py（matplotlib）生成
#   Figure 2  轴 × 免疫特征热图（a 通路级 / b 单基因级；标注基因重叠与轴内成员）
#   Figure 3  按营养轴中位数二分的 13 个免疫特征箱线图
#   Figure 4  营养轴效果量森林图（a Spearman rho + 95% bootstrap CI；
#                                b 高/低组中位数差 + 95% bootstrap CI）
#   Figure 5  （a）100 次随机分割的留出集 C-index 分布与入模特征构成
#             （b）单次 7:3 分割留出集的 KM 曲线
#
# 【2026-10-16 A6-0.2/0.3/0.4 版式合规化】依据 Discover Oncology 官方投稿指南
#   （journal 12672，条款摘录见 work/do_guidelines_excerpt.txt）：
#   - 图内一律不放标题与图注（"Do not include titles or captions within your
#     illustrations"）；描述性 title/subtitle 已删除，信息移入 03_figure_legends.md；
#     仅保留 Springer 要求的面板字母 a/b/c（粗体）
#   - 字号 8–12 pt：正文文字最小 8.0 pt（ggplot size 2.85 ≈ 8.1 pt），面板字母 12 pt
#   - 画布 ≤174 mm 宽（6.85 in）、≤234 mm 高（9.2 in）；PNG 600 dpi（Combination Art 下限）
#   - 矢量版用 cairo_pdf（字体嵌入）；原 8.8 pt subtitle 已全部删除
#   - 只改输出与样式，不改任何计算逻辑与数据
#
# 只读输入，不修改任何既有结果；新增结果写入 results/figure_data/。
# 说明：Figure 4 的 95% CI 为既有新增统计量（百分位 bootstrap，1,000 次重抽样），
#       seed 固定以保证可复现。图注中已说明 CI 的计算方式。
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(data.table)
})

FIG <- PATHA$figures
FDIR <- file.path(PATHA$results, "figure_data")
dir.create(FDIR, showWarnings = FALSE, recursive = TRUE)

SEED <- PARS$seed
BOOT <- 1000

AXCOL <- c(nutrition_liver = "#B85C1E", neutrophil = "#2C6E9B", lymphocyte = "#4A7C4E")
GRPCOL <- c(Low = "#9FB8CC", High = "#C0392B")
RHO_LOW <- "#2C6E9B"; RHO_HI <- "#C0392B"

short <- function(x) {
  x <- sub("^GOBP_", "", x)
  x <- sub("^HALLMARK_", "", x)
  gsub("_", " ", x)
}

# ---------------------------------------------------------------- 数据读取
d    <- readRDS(file.path(PATHA$work, "pa_data.rds"))
sc3  <- readRDS(file.path(PATHA$work, "pa_scores.rds"))
imm  <- readRDS(file.path(PATHA$work, "pa_scores_immune.rds"))
ex   <- d$expr_log
axes <- names(d$genesets)

r1 <- fread(file.path(PATHA$results, "immune_axis_correlation.tsv"))
r2 <- fread(file.path(PATHA$results, "checkpoint_correlation.tsv"))
t5 <- fread(file.path(PATHA$results, "table5_nutrition_axis_raw.tsv"))

S  <- as.data.frame(sc3[, axes, drop = FALSE]); rownames(S) <- sc3$sample
I  <- as.data.frame(imm)
common <- intersect(rownames(S), rownames(I))
S  <- S[common, , drop = FALSE]
I  <- I[common, , drop = FALSE]
stopifnot(nrow(S) == 1011L)

feat_levels <- t5[order(diff), immune_feature]      # 由小到大 → 画图时由下到上

# ============================================================ Figure 2a
ord_a <- r1[axis == "nutrition_liver"][order(rho), immune_feature]
A <- copy(r1)
A[, feat := factor(immune_feature, levels = ord_a)]
A[, ax := factor(axis, levels = c("nutrition_liver", "neutrophil", "lymphocyte"))]
A[, lab := sprintf("%.2f", rho)]
A[, shared  := overlap_genes > 0L]   # 有基因重叠 = 定义性关联，用黑框标注

p2a <- ggplot(A, aes(ax, feat)) +
  geom_tile(aes(fill = rho), colour = "grey85", linewidth = 0.4) +
  geom_tile(data = A[shared == TRUE], aes(ax, feat),
            fill = NA, colour = "#1A1A1A", linewidth = 1.1) +
  geom_text(aes(label = lab), size = 3.1, colour = "grey10") +
  scale_fill_gradient2(low = RHO_LOW, mid = "white", high = RHO_HI,
                       midpoint = 0, limits = c(0, 0.95),
                       name = expression(Spearman~rho)) +
  scale_x_discrete(labels = c("Nutrition", "Neutrophil", "Lymphocyte")) +
  scale_y_discrete(labels = function(x) short(as.character(x))) +
  labs(title = "a") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 0),
        legend.position = "right",
        legend.key.height = unit(1.1, "cm"))

# ============================================================ Figure 2b
CK <- c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT","IDO1","SIGLEC15")
EF <- c("IFNG","GZMA","GZMB","PRF1","CXCL9","CXCL10","CXCL11","HLA-DRA","STAT1","CCR5")
pan <- c(setNames(rep("Immune checkpoint", length(CK)), CK),
         setNames(rep("Effector / inflammation", length(EF)), EF))
B <- copy(r2)
B[, panel := pan[gene]]
B[, gene2 := factor(gene, levels = B[axis == "lymphocyte"][order(rho), gene])]
B[, ax := factor(axis, levels = c("nutrition_liver", "neutrophil", "lymphocyte"))]
B[, lab := sprintf("%.2f", rho)]
B[, in_axis := mapply(function(a, g) g %in% d$genesets[[a]], as.character(axis), gene)]
B[, sig := p < 0.05]

p2b <- ggplot(B, aes(ax, gene2)) +
  geom_tile(aes(fill = rho), colour = "grey85", linewidth = 0.4) +
  geom_tile(data = B[sig == FALSE], aes(ax, gene2),
            fill = NA, colour = "grey45", linetype = "dashed", linewidth = 0.7) +
  geom_tile(data = B[in_axis == TRUE], aes(ax, gene2),
            fill = NA, colour = "#1A1A1A", linewidth = 1.1) +
  geom_text(aes(label = lab, colour = abs(rho) > 0.6), size = 2.9) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10"), guide = "none") +
  scale_fill_gradient2(low = RHO_LOW, mid = "white", high = RHO_HI,
                       midpoint = 0, limits = c(-0.1, 0.95),
                       name = expression(Spearman~rho)) +
  scale_x_discrete(labels = c("Nutrition", "Neutrophil", "Lymphocyte")) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
  labs(title = "b") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        axis.title = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(size = 8.5),
        legend.key.height = unit(1.1, "cm"))

fig2 <- p2a / p2b + plot_layout(heights = c(13, 19))
ggsave(file.path(FIG, "fig2_axis_immune_heatmap.png"), fig2,
       width = 6.8, height = 9.0, dpi = 600)
ggsave(file.path(FIG, "fig2_axis_immune_heatmap.pdf"), fig2,
       width = 6.8, height = 9.0, device = cairo_pdf)
message("[10] Figure 2 完成")

# ============================================================ Figure 3
long <- rbindlist(lapply(feat_levels, function(f) {
  data.table(feature = f, score = I[[f]], sample = rownames(I))
}))
grp <- data.table(sample = sc3$sample,
                  nutrition_liver = sc3$nutrition_liver)
grp[, grp := ifelse(nutrition_liver > median(nutrition_liver), "High", "Low")]
long <- merge(long, grp[, .(sample, grp)], by = "sample")
long[, feature := factor(feature, levels = feat_levels)]

fdr_lab <- t5[, .(feature = immune_feature,
                  fdr_txt = ifelse(fdr < 1e-300, "FDR < 1e-300",
                                   sprintf("FDR = %.1e", fdr)))]
fdr_lab[, feature := factor(feature, levels = feat_levels)]
ymax <- long[, .(top = max(score)), by = feature]
fdr_lab <- merge(fdr_lab, ymax, by = "feature")

p3 <- ggplot(long, aes(score, feature, fill = grp)) +
  geom_boxplot(outlier.size = 0.18, outlier.alpha = 0.25,
               linewidth = 0.32, width = 0.68, position = position_dodge(width = 0.72)) +
  geom_text(data = fdr_lab, aes(x = top, y = feature, label = fdr_txt),
            inherit.aes = FALSE, hjust = -0.06, size = 2.85, colour = "grey20") +
  scale_fill_manual(values = GRPCOL, name = "Nutrition axis",
                    breaks = c("Low", "High")) +
  scale_y_discrete(labels = function(x) short(as.character(x))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.30))) +
  labs(x = "ssGSEA activity score", y = NULL) +
  theme_bw(base_size = 10.5) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "bottom")

ggsave(file.path(FIG, "fig3_nutrition_group_boxplot.png"), p3,
       width = 6.8, height = 5.4, dpi = 600)
ggsave(file.path(FIG, "fig3_nutrition_group_boxplot.pdf"), p3,
       width = 6.8, height = 5.4, device = cairo_pdf)
message("[10] Figure 3 完成")

# ============================================================ Figure 4（含 bootstrap CI）
set.seed(SEED)
nobs <- nrow(S)
b_idx <- lapply(seq_len(BOOT), function(i) sample.int(nobs, nobs, replace = TRUE))

# 4a：营养轴 × 13 特征 的 rho  bootstrap
rho_boot <- vapply(b_idx, function(k) {
  xr <- rank(S$nutrition_liver[k])
  vapply(feat_levels, function(f) {
    stats::cor(xr, rank(I[[f]][k]), method = "pearson")
  }, 0)
}, numeric(length(feat_levels)))
rho_ci <- data.table(feature = feat_levels,
                     rho = r1[axis == "nutrition_liver"][match(feat_levels, immune_feature), rho],
                     lo = apply(rho_boot, 1, stats::quantile, .025, names = FALSE),
                     hi = apply(rho_boot, 1, stats::quantile, .975, names = FALSE))

# 4b：高/低组中位数差 bootstrap
g0 <- grp$grp[match(rownames(I), grp$sample)]
diff_boot <- vapply(b_idx, function(k) {
  gg <- g0[k]
  vapply(feat_levels, function(f) {
    v <- I[[f]][k]
    median(v[gg == "High"]) - median(v[gg == "Low"])
  }, 0)
}, numeric(length(feat_levels)))
diff_ci <- data.table(feature = feat_levels,
                      diff = t5[match(feat_levels, immune_feature), diff],
                      lo = apply(diff_boot, 1, stats::quantile, .025, names = FALSE),
                      hi = apply(diff_boot, 1, stats::quantile, .975, names = FALSE))

fwrite(rho_ci,  file.path(FDIR, "fig4a_nutrition_rho_ci.tsv"),  sep = "\t")
fwrite(diff_ci, file.path(FDIR, "fig4b_nutrition_diff_ci.tsv"), sep = "\t")

forest <- function(dt, xcol, xlo, xhi, xlab, ttl, vline, xlim) {
  dt <- copy(dt)
  dt[, feature := factor(feature, levels = feat_levels)]
  ggplot(dt, aes(.data[[xcol]], feature)) +
    geom_vline(xintercept = vline, linetype = 2, colour = "grey45") +
    geom_errorbarh(aes(xmin = .data[[xlo]], xmax = .data[[xhi]]),
                   height = 0.22, colour = AXCOL[["nutrition_liver"]], linewidth = 0.5) +
    geom_point(size = 1.9, colour = AXCOL[["nutrition_liver"]]) +
    scale_y_discrete(labels = function(x) short(as.character(x))) +
    coord_cartesian(xlim = xlim) +
    labs(title = ttl, x = xlab, y = NULL) +
    theme_bw(base_size = 10.5) +
    theme(panel.grid.major.y = element_blank(),
          plot.title = element_text(face = "bold", size = 12))
}

p4a <- forest(rho_ci, "rho", "lo", "hi",
              expression("Spearman rho (95% bootstrap CI)"),
              "a", 0, c(0, 0.34))
p4b <- forest(diff_ci, "diff", "lo", "hi",
              "Median difference (high minus low)",
              "b", 0, c(0, 0.075))

fig4 <- p4a / p4b
ggsave(file.path(FIG, "fig4_nutrition_effect_forest.png"), fig4,
       width = 6.8, height = 7.6, dpi = 600)
ggsave(file.path(FIG, "fig4_nutrition_effect_forest.pdf"), fig4,
       width = 6.8, height = 7.6, device = cairo_pdf)
message("[10] Figure 4 完成（bootstrap CI 已写入 results/figure_data/）")

# ============================================================ Figure 5
rr  <- fread(file.path(PATHA$results, "internal_repeated_cv.tsv"))
aud <- fread(file.path(PATHA$results, "_audit_internal_selected.tsv"))
i1  <- fread(file.path(PATHA$results, "internal_split_summary.tsv"))

p5a <- ggplot(rr, aes(c_test)) +
  geom_histogram(bins = 18, fill = "grey72", colour = "white", linewidth = 0.25) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey30", linewidth = 0.6) +
  geom_vline(xintercept = median(rr$c_test), colour = "#C0392B", linewidth = 0.8) +
  annotate("text", x = median(rr$c_test) + 0.012, y = Inf, vjust = 1.6, hjust = 0,
           label = sprintf("median = %.3f", median(rr$c_test)),
           colour = "#C0392B", size = 3.1) +
  annotate("text", x = 0.5 - 0.012, y = Inf, vjust = 1.6, hjust = 1,
           label = "chance (0.5)", colour = "grey30", size = 3.0) +
  scale_x_continuous(limits = c(0.40, 0.60)) +
  labs(title = "a", x = "Held-out C-index", y = "Number of splits") +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", size = 12))

sel <- aud[, .N, by = selected]
sel[, selected := factor(selected, levels = c("lymphocyte", "nutrition_liver", "neutrophil"))]
p5b <- ggplot(sel, aes(selected, N)) +
  geom_col(fill = "grey62", width = 0.6) +
  geom_text(aes(label = N), vjust = -0.35, size = 3.1) +
  scale_x_discrete(labels = function(x) short(x)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(title = "b", x = NULL, y = "Number of splits") +
  theme_bw(base_size = 10.5) +
  theme(panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", size = 12))

# 5c KM：复现 R/08 的单次分割（与 results/internal_split_summary.tsv 同源）
sv <- merge(copy(d$surv_luad), sc3, by = "sample")
set.seed(PARS$seed)
idx <- sample.int(nrow(sv), floor(0.7 * nrow(sv)))
tr <- sv[idx]; te <- sv[-idx]
mdl <- readRDS(file.path(PATHA$work, "pa_model.rds"))
mu <- vapply(tr[, ..axes], mean, 0); sdv <- vapply(tr[, ..axes], stats::sd, 0)
ztr <- as.data.frame(scale(tr[, ..axes], center = mu, scale = sdv))
zte <- as.data.frame(scale(te[, ..axes], center = mu, scale = sdv))
pv <- vapply(axes, function(a)
  summary(survival::coxph(survival::Surv(tr$time, tr$status) ~ ztr[[a]]))$coefficients[1, 5], 0)
keep <- names(pv)[pv <= 0.10]; if (!length(keep)) keep <- names(pv)[which.min(pv)]
if (length(keep) == 1L) {
  cf <- setNames(unname(stats::coef(survival::coxph(survival::Surv(tr$time, tr$status) ~ ztr[[keep]]))), keep)
} else {
  fit <- suppressWarnings(fit_lasso_cox(ztr[, keep, drop = FALSE], tr$time, tr$status,
                                        nfolds = PARS$lasso_nfolds, seed = PARS$seed))
  cf <- fit$coef[fit$coef != 0]
}
rtr <- risk_score(ztr[, names(cf), drop = FALSE], cf)
rte <- risk_score(zte[, names(cf), drop = FALSE], cf)
cut <- median(rtr)
te2 <- data.table(te[, .(time, status)]); te2[, grp := ifelse(rte > cut, "High", "Low")]
f <- survival::survfit(survival::Surv(time, status) ~ grp, data = te2)
sf <- summary(f, times = sort(unique(te2$time)))
kd <- data.table(time = sf$time, surv = sf$surv,
                 grp = sub("^grp=", "", as.character(sf$strata)))
p5c <- ggplot(kd, aes(time, surv, colour = grp)) +
  geom_step(linewidth = 0.75) +
  scale_colour_manual(values = GRPCOL, name = "Risk") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "c", x = "Time (days)", y = "Overall survival") +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "bottom")

fig5 <- (p5a | p5c) / p5b + plot_layout(heights = c(2, 1))
ggsave(file.path(FIG, "fig5_internal_validation.png"), fig5,
       width = 6.8, height = 7.0, dpi = 600)
ggsave(file.path(FIG, "fig5_internal_validation.pdf"), fig5,
       width = 6.8, height = 7.0, device = cairo_pdf)
message("[10] Figure 5 完成")

# ---------------------------------------------------------------- 图注素材落盘
fwrite(rho_ci,  file.path(FDIR, "fig4a_nutrition_rho_ci.tsv"), sep = "\t")
cat("\n[10] 全部图已生成：\n")
print(list.files(FIG, pattern = "^fig[2-5]"))

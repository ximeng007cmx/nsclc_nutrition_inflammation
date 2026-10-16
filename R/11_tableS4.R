# =====================================================================
# 11_tableS4.R — 导出 Table S4（a–d）特征清单
# 目的：把 Methods 正文里为避免超长而外移的基因/通路清单落成独立附表。
# 原则：本脚本**只读取已落盘的源定义与结果**，不做任何新计算、不改数值；
#       输出的每个字段都能追到 R/00_setup_pathA.R、R/06、R/07 或 results/*.tsv。
# 输出（results/tables/）：
#   TableS4a_axis_genes.csv           3 条先验轴 × 16 基因
#   TableS4b_immune_features.csv      13 条 MSigDB 免疫特征（含与三轴的基因重叠数）
#   TableS4c_checkpoint_effector.csv  19 个免疫检查点 / 效应转录本
#   TableS4d_msigdb6_pathways.csv     S1 方案的 6 条 MSigDB 通路
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

TDIR <- file.path(PATHA$results, "tables")
dir.create(TDIR, showWarnings = FALSE, recursive = TRUE)

d  <- readRDS(file.path(PATHA$work, "pa_data.rds"))
ex <- d$expr_log
SYM <- rownames(ex)

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

label <- function(id) {
  s <- sub("^(GOBP|HALLMARK)_", "", id)
  s <- gsub("_", " ", s)
  s <- tolower(s)
  # 细胞类型/因子的首字母缩写恢复大写，避免 "Regulatory t cell" 这类写法
  s <- gsub("\\bt cell\\b", "T cell", s)
  s <- gsub("\\bb cell\\b", "B cell", s)
  s <- gsub("\\bnk cell\\b", "NK cell", s)
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}
coll_of <- function(id) if (grepl("^HALLMARK_", id)) "Hallmark" else "GO:BP"

# ---------------- S4-a 三条先验轴基因清单 ----------------
a4a <- rbindlist(lapply(names(GENES), function(ax) {
  g <- GENES[[ax]]
  data.table(axis = ax, gene_symbol = g,
             detected_in_expression_matrix = ifelse(g %in% SYM, "Y", "N"))
}))
a4a[, order_in_set := seq_len(.N), by = axis]
setcolorder(a4a, c("axis", "order_in_set", "gene_symbol",
                   "detected_in_expression_matrix"))
fwrite(a4a, file.path(TDIR, "TableS4a_axis_genes.csv"))
message(sprintf("[S4a] 3 轴 × %d 基因 = %d 行；未检出 %d 个",
                unique(a4a[, .N, by = axis]$N)[1], nrow(a4a),
                sum(a4a$detected_in_expression_matrix == "N")))

# ---------------- S4-b 13 条 MSigDB 免疫特征 ----------------
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
raw_sets <- c(read_gmt(GOBP, IMM_GO), read_gmt(HALL, IMM_HALL))
imm_sets <- lapply(raw_sets, function(g) intersect(g, SYM))

# 与三轴的基因重叠数：直接复用 results/immune_axis_correlation.tsv 的 overlap_genes
ov <- fread(file.path(PATHA$results, "immune_axis_correlation.tsv"))
ovw <- dcast(ov, immune_feature ~ axis, value.var = "overlap_genes")

a4b <- data.table(
  immune_feature  = names(imm_sets),
  collection      = vapply(names(imm_sets), coll_of, ""),
  feature_label   = vapply(names(imm_sets), label, ""),
  set_size_in_gmt = vapply(names(imm_sets), function(k) length(raw_sets[[k]]), 0L),
  set_size_in_matrix = vapply(imm_sets, length, 0L))
a4b <- merge(a4b, ovw, by = "immune_feature", all.x = TRUE, sort = FALSE)
setnames(a4b, c("nutrition_liver", "neutrophil", "lymphocyte"),
         c("overlap_nutrition_liver", "overlap_neutrophil", "overlap_lymphocyte"))
a4b[, shared_with_any_axis := ifelse(
  overlap_nutrition_liver + overlap_neutrophil + overlap_lymphocyte > 0, "Y", "N")]
setcolorder(a4b, c("immune_feature", "feature_label", "collection",
                   "set_size_in_gmt", "set_size_in_matrix",
                   "overlap_nutrition_liver", "overlap_neutrophil",
                   "overlap_lymphocyte", "shared_with_any_axis"))
setorder(a4b, collection, -set_size_in_matrix)
fwrite(a4b, file.path(TDIR, "TableS4b_immune_features.csv"))
message(sprintf("[S4b] %d 条免疫特征；规模 %d–%d；与任意轴共享基因者 %d 条",
                nrow(a4b), min(a4b$set_size_in_matrix), max(a4b$set_size_in_matrix),
                sum(a4b$shared_with_any_axis == "Y")))

# ---------------- S4-c 19 个免疫检查点 / 效应转录本 ----------------
CHECKPOINT <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2",
                "TIGIT", "IDO1", "SIGLEC15")
EFFECTOR   <- c("IFNG", "GZMA", "GZMB", "PRF1", "CXCL9", "CXCL10", "CXCL11",
                "HLA-DRA", "STAT1", "CCR5")
a4c <- rbindlist(list(
  data.table(panel = "Immune checkpoint", order_in_panel = seq_along(CHECKPOINT),
             gene_symbol = CHECKPOINT),
  data.table(panel = "Immune effector", order_in_panel = seq_along(EFFECTOR),
             gene_symbol = EFFECTOR)))
a4c[, detected_in_expression_matrix := ifelse(gene_symbol %in% SYM, "Y", "N")]
fwrite(a4c, file.path(TDIR, "TableS4c_checkpoint_effector.csv"))
message(sprintf("[S4c] 检查点 %d + 效应 %d = %d 行；全部检出 = %s",
                length(CHECKPOINT), length(EFFECTOR), nrow(a4c),
                all(a4c$detected_in_expression_matrix == "Y")))

# ---------------- S4-d S1 方案的 6 条 MSigDB 通路 ----------------
S1_HALL <- c("HALLMARK_COAGULATION", "HALLMARK_INFLAMMATORY_RESPONSE",
             "HALLMARK_ALLOGRAFT_REJECTION")
S1_GOBP <- c("GOBP_ACUTE_PHASE_RESPONSE", "GOBP_NEUTROPHIL_DEGRANULATION",
             "GOBP_ADAPTIVE_IMMUNE_RESPONSE")
s1raw <- c(read_gmt(HALL, S1_HALL), read_gmt(GOBP, S1_GOBP))
a4d <- data.table(
  pathway = names(s1raw),
  collection = vapply(names(s1raw), coll_of, ""),
  pathway_label = vapply(names(s1raw), label, ""),
  set_size_in_gmt = vapply(s1raw, length, 0L),
  set_size_in_matrix = vapply(s1raw, function(g) length(intersect(g, SYM)), 0L))
setorder(a4d, set_size_in_matrix)
fwrite(a4d, file.path(TDIR, "TableS4d_msigdb6_pathways.csv"))
message(sprintf("[S4d] %d 条通路；交集后规模 %d–%d",
                nrow(a4d), min(a4d$set_size_in_matrix), max(a4d$set_size_in_matrix)))

message("\n[11] Table S4（a–d）已写入 ", TDIR)

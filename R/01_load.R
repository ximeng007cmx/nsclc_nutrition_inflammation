# =====================================================================
# 01_load.R — 读取表达矩阵 + 元数据 + 临床表，构造 LUAD(训练) / LUSC(验证)
# 输出 work/pa_data.rds：
#   expr_log  : log2(FPKM+1) 矩阵（symbol 级，仅 -01 原发肿瘤样本）
#   surv_luad : LUAD 样本级生存表（样本条码 + 临床）
#   surv_lusc : LUSC 同上
# =====================================================================
source("R/00_setup_pathA.R")
source("R/pathA_core.R")

message("[01] 读取表达矩阵 …")
expr_raw <- read_expr(PATHA$expr)
message(sprintf("     %d 基因 × %d 样本（raw, ENSG 级）",
                nrow(expr_raw), ncol(expr_raw)))

meta <- fread(PATHA$expr_meta)
# 仅保留原发肿瘤（-01）
meta_t <- meta[grepl("-01$", sample)]
expr_t <- expr_raw[, colnames(expr_raw) %in% meta_t$sample, drop = FALSE]
message(sprintf("     原发肿瘤 %d 样本（LUAD+LUSC）", ncol(expr_t)))

message("[01] ENSG -> symbol …")
pm <- read_probemap(PATHA$probemap)
expr_sym <- map_to_symbol(expr_t, pm)
expr_log <- toil_fix(expr_sym)
message(sprintf("     %d 基因（symbol 级）", nrow(expr_log)))

# 基因集覆盖检查（真实记录缺失，不硬造）
genes_avail <- rownames(expr_log)
gs <- lapply(GENES, function(g) intersect(g, genes_avail))
gs_missing <- lapply(GENES, function(g) setdiff(g, genes_avail))
for (nm in names(GENES)) {
  message(sprintf("     基因集 %-15s 命中 %2d/%2d%s", nm, length(gs[[nm]]),
                  length(GENES[[nm]]),
                  if (length(gs_missing[[nm]]))
                    paste0("，缺失: ", paste(gs_missing[[nm]], collapse = ","))
                  else ""))
}

message("[01] 临床表合并 …")
cl <- fread(PATHA$clinical)
surv_luad <- build_surv_table(meta[project == "TCGA-LUAD"], cl, PARS$endpoint)
surv_lusc <- build_surv_table(meta[project == "TCGA-LUSC"], cl, PARS$endpoint)
surv_luad[, stage_late := stage_binary(stage)]
surv_lusc[, stage_late := stage_binary(stage)]
surv_luad[, gender := factor(gender)]
surv_lusc[, gender := factor(gender)]

message(sprintf("     LUAD 生存表 %d 例（死亡 %d）",
                nrow(surv_luad), sum(surv_luad$status == 1)))
message(sprintf("     LUSC 生存表 %d 例（死亡 %d）",
                nrow(surv_lusc), sum(surv_lusc$status == 1)))

saveRDS(list(expr_log = expr_log,
             genesets = gs, genesets_missing = gs_missing,
             surv_luad = surv_luad, surv_lusc = surv_lusc),
        file.path(PATHA$work, "pa_data.rds"))
message("[01] 完成 -> work/pa_data.rds")

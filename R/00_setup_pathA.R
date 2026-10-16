# =====================================================================
# 00_setup_pathA.R — 路径 A（营养-炎症分子特征 → NSCLC 预后模型）全局配置
# 数据真实性约定：
#   * 表达数据 = UCSC Xena Toil hub 的 TCGA pan-cancer RSEM FPKM（真实下载）
#   * 临床数据 = TCGA-CDR (Liu J et al., Cell 2018; PMID 29625055)
#   * 基因集为「ALI 三组分的已知标志基因」代理（见 GENES 注释），
#     属于先验定义，不是数据驱动筛选；正文表述必须用
#     「营养-炎症分子特征」，不得说「等价于 ALI」。
# =====================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ---- 路径 ----
PATHA <- list(
  dir      = "D:/WorkBuddy_Y02/out_ali_directions/pathA",
  expr     = "D:/WorkBuddy_Y02/out_ali_directions/pathA/work/expr_luad_lusc.tsv",
  expr_meta= "D:/WorkBuddy_Y02/out_ali_directions/pathA/work/expr_luad_lusc.tsv.meta.tsv",
  probemap = "D:/WorkBuddy_Y02/out_ali_directions/_gencode_probemap.gz",  # 实为纯文本
  clinical = "D:/WorkBuddy_Y02/out_ali_directions/pathA/data/tcga_cdr_clinical.csv",
  work     = "D:/WorkBuddy_Y02/out_ali_directions/pathA/work",
  results  = "D:/WorkBuddy_Y02/out_ali_directions/pathA/results",
  figures  = "D:/WorkBuddy_Y02/out_ali_directions/pathA/figures",
  logs     = "D:/WorkBuddy_Y02/out_ali_directions/pathA/logs"
)
for (d in c(PATHA$work, PATHA$results, PATHA$figures, PATHA$logs))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- 分析参数 ----
PARS <- list(
  endpoint   = "OS",          # 主终点 OS；可切 DSS / PFI
  train_proj = "TCGA-LUAD",   # 训练集（原发肿瘤 -01）
  val_proj   = "TCGA-LUSC",   # 外部泛化（原发肿瘤 -01）
  lasso_nfolds = 10,
  lasso_alpha  = 1,
  boot_n       = 500,         # C-index bootstrap CI
  seed         = 20260101
)

# ---- 营养-炎症轴基因集（先验定义；符号级；缺失基因自动剔除并记录） ----
GENES <- list(
  # 轴 1 肝源性营养轴（白蛋白/转运蛋白类，对应 ALI 的白蛋白组分）
  nutrition_liver = c("ALB","TF","RBP4","APOA1","APOA2","APOB","GC","TTR",
                      "HPX","SERPINA1","AHSG","ORM1","FGA","FGG","ITIH2","PLG"),
  # 轴 2 中性粒细胞炎症轴（对应 ALI 的中性粒细胞计数组分）
  neutrophil = c("ELANE","MPO","LTF","CEACAM8","S100A8","S100A9","S100A12",
                 "MMP9","PADI4","CXCR2","FPR1","DEFA4","CAMP","LCN2","ARG1","OLFM4"),
  # 轴 3 淋巴细胞/细胞毒 T 轴（对应 ALI 的淋巴细胞计数组分）
  lymphocyte = c("CD3D","CD3E","CD3G","CD8A","CD8B","LCK","GZMB","GZMK",
                 "PRF1","IL7R","CCL5","CXCL10","CXCR6","CD27","CXCL13","HLA-DRA")
)

message("[00_setup] 配置完成。基因集规模: ",
        paste(sprintf("%s=%d", names(GENES), vapply(GENES, length, 0L)), collapse=", "))

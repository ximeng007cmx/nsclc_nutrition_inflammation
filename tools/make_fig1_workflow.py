# -*- coding: utf-8 -*-
"""
tools/make_fig1_workflow.py — Figure 1：研究设计与分析流程图

【2026-10-16 A6-0.2 版式重排】
按 Discover Oncology 官方投稿指南（journal 12672，见 work/do_guidelines_excerpt.txt）重排：
  - 字体：sans（8-12 pt）；标题 10 pt、盒首行 8.6 pt、盒正文 8.0 pt，全部落在区间内
  - 画布宽 6.85 in（= 174 mm，Springer 单栏文本区宽度），高 7.5 in（= 190.5 mm < 234 mm）
  - 旧版 15.2 in 宽的横向布局缩到 174 mm 后正文仅 ~4 pt，不合规，故重排为四段纵向布局
  - 无图内标题/图注/脚注（原底部脚注句已由 Methods 2.3/2.5 正文覆盖，图注文件另录）
输出：figures/fig1_workflow.png（600 dpi）、figures/fig1_workflow.pdf（字体嵌入）

数字来源核对（与 A2/A3 一致，未改动任何数字）：
  10,535 样本列      -> gzip 读表头实测（tools/audit_methods_facts.R 同源）
  1,011 原发瘤       -> work/pa_data.rds 实测
  991 / 500 / 491    -> results/tables/Table1_baseline.csv
  3 x 16 基因        -> R/00_setup_pathA.R GENES
  13 免疫特征 / 19 转录本 -> results/immune_axis_correlation.tsv、checkpoint_correlation.tsv
  C-index 0.494 / 0.542 -> results/internal_split_summary.tsv、internal_repeated_cv.tsv（2026-10-18 重跑后值）
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams["font.family"] = "DejaVu Sans"   # sans-serif；投稿包可另存 Arial/Helvetica 版

OUT_PNG = r"D:/WorkBuddy_Y02/out_ali_directions/pathA/figures/fig1_workflow.png"
OUT_PDF = r"D:/WorkBuddy_Y02/out_ali_directions/pathA/figures/fig1_workflow.pdf"

# 6.85 in = 174 mm（Springer 单栏文本区宽度）；7.5 in = 190.5 mm < 234 mm 上限
fig, ax = plt.subplots(figsize=(6.85, 8.25))   # 174 x 210 mm
ax.set_xlim(0, 180)
ax.set_ylim(0, 220)
ax.axis("off")

C_DATA = "#EDF2F7"; E_DATA = "#7A8CA3"
C_AX = "#FDF0E4";  E_AX = "#B85C1E"
C_IMM = "#E8F1F8"; E_IMM = "#2C6E9B"
C_MOD = "#EAF3EC"; E_MOD = "#4A7C4E"
C_RES = "#FFF8E1"; E_RES = "#B08A2A"

FS_TITLE = 10.0   # 栏目标题
FS_HEAD = 8.6     # 盒内标题行（粗体，可多行）
FS_BODY = 8.0     # 盒内正文（>= 8 pt 下限）


TH = 4.0   # 标题行高（数据坐标单位，实测校准）
BH = 5.5   # 正文行高（数据坐标单位，实测校准）
GAP = 2.0  # 标题与正文间距


def box(x, y, w, h, title, body, fc, ec, lw=1.2):
    """title: 粗体标题，可含换行；body: 正文，可含换行（可为空串）。内容在盒内垂直居中。"""
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.5,rounding_size=1.4",
                                linewidth=lw, edgecolor=ec, facecolor=fc, zorder=2))
    nt = title.count("\n") + 1
    nb = body.count("\n") + 1 if body else 0
    total = TH * nt + (GAP + BH * nb if nb else 0)
    top = y + (h + total) / 2.0
    ax.text(x + w / 2, top, title, ha="center", va="top",
            fontsize=FS_HEAD, fontweight="bold", color="#1A1A1A",
            linespacing=1.25, zorder=3)
    if body:
        ax.text(x + w / 2, top - TH * nt - GAP, body, ha="center", va="top",
                fontsize=FS_BODY, color="#333333", linespacing=1.4, zorder=3)


def arrow(x1, y1, x2, y2, color="#5A6B7B", lw=1.3, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=11, linewidth=lw, color=color,
                                 zorder=1, shrinkA=1, shrinkB=1))


def band_label(x, y, text, color):
    ax.text(x, y, text, ha="left", va="center",
            fontsize=FS_TITLE, fontweight="bold", color=color)


COL = [(8, 56), (67, 56), (125, 47)]   # 三列 (x, w)

# ---------------- Band A：数据源（y 170-190） ----------------
band_label(8, 214, "Data sources", "#1A1A1A")
box(COL[0][0], 182, COL[0][1], 30,
    "TCGA RNA-seq",
    "Toil reprocessed\ncompendium\n10,535 sample columns",
    C_DATA, E_DATA)
box(COL[1][0], 182, COL[1][1], 30,
    "GDC barcodes",
    "primary solid tumours\n(sample-type code 01)",
    C_DATA, E_DATA)
box(COL[2][0], 182, COL[2][1], 30,
    "TCGA-CDR",
    "curated OS, DSS\nand PFI",
    C_DATA, E_DATA)

# ---------------- Band B：队列汇总（y 146-164） ----------------
box(8, 150, 164, 22,
    "Cohort assembly",
    "1,011 primary tumours\n"
    "991 with survival annotation (500 LUAD, 491 LUSC)",
    C_DATA, E_DATA)
for x, w in COL:
    arrow(x + w / 2, 182, x + w / 2, 172, color=E_DATA)

# ---------------- Band C：主目标（关联分析，y 66-124） ----------------
band_label(8, 143, "Association analysis (primary aim)", E_AX)
ax.plot([8, 172], [139, 139], color=E_AX, lw=1.8, solid_capstyle="round")

box(COL[0][0], 71, COL[0][1], 64,
    "Nutrition\u2013\ninflammation axes",
    "3 a priori gene sets\n(16 genes each)\nssGSEA activity scores",
    C_AX, E_AX)

box(COL[1][0], 71, COL[1][1], 64,
    "Immune\nmicroenvironment",
    "13 MSigDB\nimmune features\n19 checkpoint /\neffector transcripts\n"
    "gene overlap\nannotated per pair\n(0 vs \u22651 shared)",
    C_IMM, E_IMM)

box(COL[2][0], 71, COL[2][1], 64,
    "Association\nresults",
    "Spearman rho\n+ BH FDR\nmedian-split\ncomparison\n(Wilcoxon rank-sum)\n"
    "nutrition axis\ntracks myeloid\ninflammation with\nno shared genes",
    C_RES, E_RES)

arrow(90, 150, 36, 135, color=E_AX)
arrow(64, 103, 67, 103, color=E_IMM)
arrow(123, 103, 125, 103, color=E_RES)
arrow(36, 71, 36, 56, color=E_MOD)

# ---------------- Band D：次目标（预后评估，y 10-48） ----------------
band_label(8, 64.5, "Prognostic evaluation (secondary aim)", E_MOD)
ax.plot([8, 172], [60.5, 60.5], color=E_MOD, lw=1.8, solid_capstyle="round")

box(COL[0][0], 14, COL[0][1], 42,
    "LASSO\u2013Cox composite",
    "lambda.min\npre-specified zero-\ncoefficient fallback\n"
    "training-set median\ncut point",
    C_MOD, E_MOD)

box(COL[1][0], 14, COL[1][1], 42,
    "Internal validation",
    "random 7:3 split\ntraining-set-only\nscreening and\nstandardisation\n"
    "100 repeated splits",
    C_MOD, E_MOD)

box(COL[2][0], 14, COL[2][1], 42,
    "Prognostic result",
    "endpoint sensitivity\nand feature schemes\n"
    "no prognostic\ndiscrimination\n(C-index 0.494;\nmedian 0.542)",
    C_RES, E_RES)

arrow(64, 35, 67, 35, color=E_MOD)
arrow(123, 35, 125, 35, color=E_RES)

fig.savefig(OUT_PNG, dpi=600, bbox_inches="tight", facecolor="white")
fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
print("OK fig1_workflow.png (600 dpi) / .pdf")

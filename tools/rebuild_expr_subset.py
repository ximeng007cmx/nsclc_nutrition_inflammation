"""按条码名从原始 Xena 矩阵重建 work/expr_luad_lusc.tsv。

背景（2026-10-17）：旧子集的列标签与列内容整体错位 —— 用 400 基因指纹在原始
矩阵中反查，12 个随机列无一匹配自身条码（10 个以 r=1.00000 匹配到别的样本），
且原始矩阵的性染色体基因信号干净（XIST 女 5.21 / 男 -6.08，子集里两性均 ~3.0）。
因此子集不可用，只能按名重建。

做法：从 data/tcga_RSEM_gene_fpkm.gz 流式读取，按 meta 的 1,120 个条码名取列，
保持原始字符串（不做浮点round-trip），基因行顺序与原始一致。
"""
import csv
import gzip
import os
import sys
import time

RAW = "data/tcga_RSEM_gene_fpkm.gz"
META = "work/expr_luad_lusc.tsv.meta.tsv"
OUT = "work/expr_luad_lusc.tsv"

meta = []
with open(META, encoding="utf-8-sig") as f:
    h = f.readline().rstrip("\n").split("\t")
    si = h.index("sample")
    for row in csv.reader(f, delimiter="\t"):
        if row and row[si]:
            meta.append(row[si])

t0 = time.time()
with gzip.open(RAW, "rt") as fin:
    header = fin.readline().rstrip("\n").split("\t")
    raw_samples = header[1:]
    pos = {s: i for i, s in enumerate(raw_samples)}
    missing = [s for s in meta if s not in pos]
    if missing:
        sys.exit("meta 中 %d 个条码不在原始矩阵：%s" % (len(missing), missing[:3]))
    take = [pos[s] for s in meta]           # meta 顺序 -> 原始列索引

    with open(OUT, "w", newline="") as fout:
        fout.write("gene\t" + "\t".join(meta) + "\n")
        n = 0
        buf = []
        for line in fin:
            p = line.rstrip("\n").split("\t")
            buf.append("\t".join([p[0]] + [p[i + 1] for i in take]) + "\n")
            n += 1
            if len(buf) >= 2000:
                fout.writelines(buf)
                buf = []
            if n % 20000 == 0:
                print("  %d genes, %.1f min" % (n, (time.time() - t0) / 60), flush=True)
        fout.writelines(buf)

print("written %s  %.1f MB  %.1f min" % (OUT, os.path.getsize(OUT) / 1e6, (time.time() - t0) / 60))

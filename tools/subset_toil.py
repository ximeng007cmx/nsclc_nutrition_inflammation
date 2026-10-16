# -*- coding: utf-8 -*-
"""把 Toil pan-cancer 表达矩阵（gz，全癌种）流式子集为 LUAD/LUSC。

内存友好：逐行读、只保留目标样本列。
输出：
  <out>.tsv       基因 × 样本 表达矩阵（首列 sample = 基因 ENSG id，其余列 = 4段样本条码）
  <out>.meta.tsv  样本元数据（sample / project / patient，每样本一行）
用法：
  python subset_toil.py <matrix.gz> <barcodes.tsv> <out.tsv> \
      [--project TCGA-LUAD,TCGA-LUSC] [--suffix 01,11]
"""
import argparse
import gzip
import sys
import time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix")
    ap.add_argument("barcodes")
    ap.add_argument("out")
    ap.add_argument("--project", default="TCGA-LUAD,TCGA-LUSC")
    ap.add_argument("--suffix", default="01,11")
    a = ap.parse_args()

    want_proj = set(a.project.split(","))
    want_sfx = {s.zfill(2) for s in a.suffix.split(",")}

    # 1) 目标 4 段条码。full_sample_id 形如 TCGA-44-3918-01A（第4段含 vial 字母），
    #    矩阵列名形如 TCGA-44-3918-01 → 第4段取前 2 个字符。同键多 vial 取首个。
    want_all = {}
    n_proj = {}
    with open(a.barcodes, encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\r\n").split("\t")
        ix = {k: i for i, k in enumerate(hdr)}
        for ln in fh:
            p = ln.rstrip("\r\n").split("\t")
            if len(p) < 3:
                continue
            proj, pat = p[ix["project"]], p[ix["patient"]]
            fsid = p[ix["full_sample_id"]]
            if proj not in want_proj:
                continue
            segs = fsid.split("-")
            if len(segs) < 4:
                continue
            key = "%s-%s-%s-%s" % (segs[0], segs[1], segs[2], segs[3][:2])
            sfx = segs[3][:2]
            if sfx not in want_sfx:
                continue
            if key in want_all:
                continue  # 同一 4 段键去重（保留首个 vial）
            want_all[key] = (proj, pat)
            n_proj[proj] = n_proj.get(proj, 0) + 1

    print("目标样本 %d 个（%s）" % (len(want_all),
          " / ".join("%s %d" % (k, v) for k, v in sorted(n_proj.items()))), flush=True)

    # 2) 流式解压，逐行写
    t0 = time.time()
    n_gene = 0
    found = set()
    with gzip.open(a.matrix, "rt") as fin, \
         open(a.out, "w", encoding="utf-8", newline="") as fout, \
         open(a.out + ".meta.tsv", "w", encoding="utf-8", newline="") as fmeta:
        header = fin.readline().rstrip("\n").split("\t")
        cols = header[1:]
        idx = [i for i, c in enumerate(cols) if c in want_all]
        sel_ids = [cols[i] for i in idx]
        found = set(sel_ids)
        fout.write("gene\t" + "\t".join(sel_ids) + "\n")
        fmeta.write("sample\tproject\tpatient\n")
        for c in sel_ids:
            proj, pat = want_all[c]
            fmeta.write("%s\t%s\t%s\n" % (c, proj, pat))
        fmeta.flush()

        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            fout.write(parts[0] + "\t" + "\t".join(parts[i] for i in idx) + "\n")
            n_gene += 1
            if n_gene % 10000 == 0:
                print("  已处理 %d 基因… %.0fs" % (n_gene, time.time() - t0), flush=True)

    print("完成：%d 个基因 × %d 个样本，用时 %.0fs" % (n_gene, len(sel_ids), time.time() - t0))
    missing = len(want_all) - len(found)
    print("未在矩阵中找到的目标样本：%d" % missing)
    if missing:
        miss = sorted(set(want_all) - found)
        print("  示例:", miss[:5])
    return 0


if __name__ == "__main__":
    sys.exit(main())

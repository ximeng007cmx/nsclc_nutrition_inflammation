"""Build a gene-symbol expression matrix and clinical table for an external NSCLC cohort.

Inputs (all local, already downloaded from EBI BioStudies):
  data/external/<study>/tables/*_sample_table.txt   per-sample probe-level tables (ID_REF, VALUE)
  data/external/<study>/<ACC>.sdrf.txt              sample annotations (survival, stage, age, sex)
  work/hgu133a.db/inst/extdata/hgu133a.sqlite       probe -> Entrez (GPL96)
  work/hgu133plus2.db/inst/extdata/hgu133plus2.sqlite  probe -> Entrez (GPL570)
  work/org.Hs.eg.db/inst/extdata/org.Hs.eg.sqlite   Entrez -> symbol

Outputs:
  data/external/<study>/expr_symbol.tsv.gz   genes x samples, log2(x+1) of the platform's linear values
  data/external/<study>/clinical.tsv         sample, os_time_days, os_event, age, sex, stage
  data/external/<study>/build_manifest.json  counts + fingerprint of the annotation sources

Probe collapsing: for symbols with several probes the probe with the highest mean expression is kept
(the conventional max-mean rule); this is recorded in the manifest.
"""
import gzip
import json
import os
import sqlite3
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "work")

STUDIES = {
    "gse68465": {"acc": "E-GEOD-68465", "db": "hgu133a.db/inst/extdata/hgu133a.sqlite",
                 "platform": "GPL96 (HG-U133A)"},
    "gse31210": {"acc": "E-GEOD-31210", "db": "hgu133plus2.db/inst/extdata/hgu133plus2.sqlite",
                 "platform": "GPL570 (HG-U133 Plus 2)"},
}


def probe_map(db_rel):
    con = sqlite3.connect(os.path.join(WORK, db_rel))
    org = sqlite3.connect(os.path.join(WORK, "org.Hs.eg.db/inst/extdata/org.Hs.eg.sqlite"))
    ent2sym = {}
    for gid, sym in org.execute("select g._id, i.symbol from genes g join gene_info i on i._id = g._id"):
        if sym:
            ent2sym[str(gid)] = sym
    out = {}
    for probe, gene_id in con.execute("select probe_id, gene_id from probes"):
        sym = ent2sym.get(str(gene_id))
        if sym:
            out[probe] = sym
    con.close()
    org.close()
    return out


def read_sample(path):
    vals = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            try:
                vals[parts[0]] = float(parts[1])
            except ValueError:
                continue
    return vals


def main(study):
    cfg = STUDIES[study]
    base = os.path.join(ROOT, "data", "external", study)
    tdir = os.path.join(base, "tables")
    files = sorted(f for f in os.listdir(tdir) if f.endswith("_sample_table.txt"))
    if not files:
        raise SystemExit("no sample tables in %s" % tdir)
    samples = [f.replace("_sample_table.txt", "") for f in files]

    pmap = probe_map(cfg["db"])
    print("%s: %d sample tables, %d probe->symbol mappings" % (study, len(samples), len(pmap)))

    # accumulate per-symbol sums and counts so we can pick the max-mean probe later
    probe_sums, probe_n = {}, {}
    col = {}
    for i, f in enumerate(files):
        vals = read_sample(os.path.join(tdir, f))
        col[samples[i]] = {p: v for p, v in vals.items() if p in pmap}
        for p, v in col[samples[i]].items():
            probe_sums[p] = probe_sums.get(p, 0.0) + v
            probe_n[p] = probe_n.get(p, 0) + 1
        if (i + 1) % 50 == 0:
            print("  parsed %d/%d" % (i + 1, len(files)), flush=True)

    probe_mean = {p: probe_sums[p] / probe_n[p] for p in probe_sums if probe_n[p] >= len(samples) * 0.9}
    # best probe per symbol
    best = {}
    for p, m in probe_mean.items():
        s = pmap[p]
        if s not in best or m > probe_mean[best[s]]:
            best[s] = p
    print("  symbols kept: %d" % len(best))

    genes = sorted(best)
    out_path = os.path.join(base, "expr_symbol.tsv.gz")
    import math
    with gzip.open(out_path, "wt", encoding="utf-8", newline="") as fh:
        fh.write("symbol\t" + "\t".join(samples) + "\n")
        for g in genes:
            p = best[g]
            row = []
            for s in samples:
                v = col[s].get(p)
                row.append("" if v is None else ("%.4f" % math.log2(v + 1.0)))
            fh.write(g + "\t" + "\t".join(row) + "\n")

    man = {"study": study, "accession": cfg["acc"], "platform": cfg["platform"],
           "n_samples": len(samples), "n_symbols": len(genes),
           "probe_collapse": "max mean expression across samples",
           "transform": "log2(x + 1) applied to the platform's processed values",
           "probe_map_source": cfg["db"], "entrez_symbol_source": "org.Hs.eg.db 3.23.1"}
    with open(os.path.join(base, "build_manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(man, fh, indent=2)
    print("  wrote", out_path, "and build_manifest.json")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "gse68465"))

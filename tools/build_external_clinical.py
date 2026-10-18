"""Parse the SDRF of an external cohort into clinical.tsv (sample, OS, age, sex, stage).

SDRF is the ArrayExpress tab-separated sample/data relationship file; column names repeat
(Term Source REF blocks are interleaved), so columns are located by exact header text and the
first occurrence of each characteristic is used.

Unit handling: GSE68465 reports months, converted with 30.4375 days/month and the conversion is
recorded in the output header comment; GSE31210 reports days directly.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SPECS = {
    "gse68465": {
        "acc": "E-GEOD-68465",
        "time_col": "Characteristics [months_to_last_contact_or_death]",
        "time_unit": "months", "to_days": 30.4375,
        "status_col": "Characteristics [vital_status]",
        "dead_values": {"dead", "deceased", "1"},
    },
    "gse31210": {
        "acc": "E-GEOD-31210",
        "time_col": "FactorValue [DAYS BEFORE DEATH/CENSOR]",
        "time_unit": "days", "to_days": 1.0,
        "status_col": "FactorValue [DEATH]",
        "dead_values": {"dead", "1", "yes", "death"},
    },
}


def find_col(header, name):
    want = name.strip().lower()
    for i, h in enumerate(header):
        if h.strip().lower() == want:
            return i
    return None


def main(study):
    spec = SPECS[study]
    base = os.path.join(ROOT, "data", "external", study)
    sdrf = os.path.join(base, spec["acc"] + ".sdrf.txt")
    with open(sdrf, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]

    i_src = find_col(header, "Source Name")
    i_t = find_col(header, spec["time_col"])
    i_s = find_col(header, spec["status_col"])
    def best_col(cands):
        """Among candidate column names pick the one with most non-empty cells."""
        idxs = [i for i, h in enumerate(header) if h.strip().lower() in cands]
        if not idxs:
            return None
        return max(idxs, key=lambda i: sum(1 for r in rows if i < len(r) and r[i].strip()))

    i_age = best_col(("characteristics [age]", "characteristics[age (years)]",
                      "factorvalue [age]", "factorvalue [age (years)]"))
    i_sex = best_col(("characteristics [sex]", "characteristics[gender]",
                      "factorvalue [gender]", "factorvalue [sex]"))
    i_stage = (find_col(header, "Characteristics [disease_stage]")
               or find_col(header, "FactorValue [PATHOLOGICAL STAGE]"))
    i_tissue = find_col(header, "FactorValue [TISSUE]")
    if None in (i_src, i_t, i_s):
        raise SystemExit("required SDRF columns missing: %s / %s / %s" % (spec["time_col"], spec["status_col"], "Source Name"))
    print("columns: source=%s time=%s status=%s age=%s sex=%s stage=%s tissue=%s" %
          (i_src, i_t, i_s, i_age, i_sex, i_stage, i_tissue))

    out = os.path.join(base, "clinical.tsv")
    n = 0
    with open(out, "w", encoding="utf-8", newline="") as fh:
        fh.write("# os_time_days derived from %s (unit=%s, factor=%s)\n" % (spec["time_col"], spec["time_unit"], spec["to_days"]))
        fh.write("sample\tos_time_days\tos_event\tage\tsex\tstage\ttissue\n")
        for r in rows:
            def get(i):
                return r[i].strip() if i is not None and i < len(r) else ""
            t = get(i_t)
            st = get(i_s).lower()
            try:
                days = float(t) * spec["to_days"]
            except ValueError:
                continue
            if days <= 0:
                continue
            event = "1" if st in spec["dead_values"] else "0"
            fh.write("\t".join([get(i_src).split()[0] if get(i_src) else "", "%.4f" % days, event, get(i_age), get(i_sex),
                                get(i_stage), get(i_tissue)]) + "\n")
            n += 1
    print("wrote %s with %d samples" % (out, n))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "gse68465"))

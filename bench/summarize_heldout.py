"""Pool held-out benchmark results over seeds, with 95% Wilson confidence intervals.

Input is a TSV with the columns: seed, mode, tech, coverage, klass, truth, tp, calls, fp.
`tp` is the number of truth variants found. `calls` is the number of calls of that class and `fp` the
number of those that are not in the truth set. Usage:

    python3 summarize_heldout.py results.tsv
"""
import argparse
import csv
import math
from collections import defaultdict

Z = 1.959964  # 95%


def wilson(k, n, z=Z):
    """The Wilson score interval for k successes in n trials. (nan, nan) when n is 0."""
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    lo, hi = centre - half, centre + half
    return (0.0 if k == 0 else max(0.0, lo), 1.0 if k == n else min(1.0, hi))


def _ratio(num, den):
    return num / den if den else math.nan


def summarize(rows):
    """Group by (mode, tech, coverage, class), pool the counts over seeds and add the intervals."""
    groups = defaultdict(list)
    for r in rows:
        groups[(r["mode"], r["tech"], int(r["coverage"]), r["klass"])].append(r)
    table = []
    for (mode, tech, coverage, klass), rs in sorted(groups.items()):
        truth = sum(int(r["truth"]) for r in rs)
        tp = sum(int(r["tp"]) for r in rs)
        calls = sum(int(r["calls"]) for r in rs)
        fp = sum(int(r["fp"]) for r in rs)
        per_seed = [_ratio(int(r["tp"]), int(r["truth"])) for r in rs if int(r["truth"])]
        s_lo, s_hi = wilson(tp, truth)
        p_lo, p_hi = wilson(calls - fp, calls)
        table.append({
            "mode": mode, "tech": tech, "coverage": coverage, "klass": klass, "seeds": len(rs),
            "truth": truth, "tp": tp, "calls": calls, "fp": fp,
            "sensitivity": _ratio(tp, truth), "sens_lo": s_lo, "sens_hi": s_hi,
            "sens_min": min(per_seed) if per_seed else math.nan, "sens_max": max(per_seed) if per_seed else math.nan,
            "precision": _ratio(calls - fp, calls), "prec_lo": p_lo, "prec_hi": p_hi,
        })
    return table


def _fmt(x):
    return "NA" if isinstance(x, float) and math.isnan(x) else (f"{x:.3f}" if isinstance(x, float) else str(x))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Pool held-out results over seeds.")
    ap.add_argument("results")
    args = ap.parse_args(argv)
    with open(args.results) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    columns = ["mode", "tech", "coverage", "klass", "seeds", "truth", "tp", "fp", "calls",
               "sensitivity", "sens_lo", "sens_hi", "sens_min", "sens_max", "precision", "prec_lo", "prec_hi"]
    print("\t".join(columns))
    for r in summarize(rows):
        print("\t".join(_fmt(r[c]) for c in columns))


if __name__ == "__main__":
    main()

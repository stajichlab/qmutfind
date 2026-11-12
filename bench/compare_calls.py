"""Compare a call set with a truth set: sensitivity and precision per variant class.

Both sets must be normalized first (`bcftools norm -f reference.fa`). A variant is
keyed by (contig, pos, ref, alt). A call matches only when all four fields are equal.
"""
import argparse
import math

from make_truth import _open_text

CLASSES = ("SNP", "SMALL_INDEL", "LARGE_INDEL")
LARGE_INDEL_MIN = 50
# ASSUMED tolerances for matching large indels. They are not derived from data. See PLAN.md.
POS_TOL = 50  # bp between the truth and call positions
SIZE_FRACTION = 0.8  # the smaller size over the larger size


def _ratio(numerator, denominator):
    return numerator / denominator if denominator else float("nan")


def _type_and_size(key, svlens):
    _contig, _pos, ref, alt = key
    if key in svlens:  # a symbolic insertion: the size is in SVLEN
        return "INS", svlens[key]
    return ("DEL" if len(ref) > len(alt) else "INS"), abs(len(ref) - len(alt))


def _match_large(truth_keys, call_keys, truth_svlen, call_svlen, pos_tol, size_frac):
    """One-to-one tolerant matching. Returns the number of matched pairs.
    Nearest pairs are matched first. A pair needs the same contig and type, a position
    within `pos_tol`, and sizes whose ratio (small over large) is at least `size_frac`."""
    pairs = []
    for t in truth_keys:
        t_type, t_size = _type_and_size(t, truth_svlen)
        for c in call_keys:
            if c[0] != t[0]:
                continue
            c_type, c_size = _type_and_size(c, call_svlen)
            if c_type != t_type or abs(c[1] - t[1]) > pos_tol:
                continue
            if min(t_size, c_size) < size_frac * max(t_size, c_size):
                continue
            pairs.append((abs(c[1] - t[1]), abs(c_size - t_size), t, c))
    used_t, used_c = set(), set()
    for _dist, _dsize, t, c in sorted(pairs, key=lambda p: p[:2]):
        if t not in used_t and c not in used_c:
            used_t.add(t)
            used_c.add(c)
    return len(used_t)


def compare_variants(truth, calls, truth_svlen=None, call_svlen=None,
                     pos_tol=POS_TOL, size_frac=SIZE_FRACTION):
    """Return {class: {truth, tp, fn, calls, fp, sensitivity, precision}}.
    `truth` and `calls` map a variant key to its class. SNPs and small indels need an
    exact match. Large indels match with a tolerance (see `_match_large`).
    `*_svlen` map a symbolic key to its SVLEN."""
    truth_svlen, call_svlen = truth_svlen or {}, call_svlen or {}
    stats = {}
    for klass in CLASSES:
        truth_keys = {k for k, c in truth.items() if c == klass}
        call_keys = {k for k, c in calls.items() if c == klass}
        if klass == "LARGE_INDEL":
            tp = _match_large(truth_keys, call_keys, truth_svlen, call_svlen, pos_tol, size_frac)
            fp = len(call_keys) - tp
        else:
            tp = len(truth_keys & set(calls))
            fp = len(call_keys - set(truth))
        stats[klass] = {
            "truth": len(truth_keys),
            "tp": tp,
            "fn": len(truth_keys) - tp,
            "calls": len(call_keys),
            "fp": fp,
            "sensitivity": _ratio(tp, len(truth_keys)),
            "precision": _ratio(len(call_keys) - fp, len(call_keys)),
        }
    return stats


def load_svlens(path):
    """Read SVLEN from INFO for records that have it: {(contig, pos, ref, alt): svlen}."""
    svlens = {}
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            if "SVLEN" in info:
                svlens[(f[0], int(f[1]), f[3], f[4])] = abs(int(info["SVLEN"].split(",")[0]))
    return svlens


def _class_from_alleles(ref, alt, large_min):
    if len(ref) == len(alt):
        return "SNP"
    return "LARGE_INDEL" if abs(len(ref) - len(alt)) >= large_min else "SMALL_INDEL"


def load_vcf(path, large_min=LARGE_INDEL_MIN):
    """Read a VCF (plain, .gz or .zst) into {(contig, pos, ref, alt): class}.
    The class comes from INFO CLASS when present, else from the allele lengths."""
    variants = {}
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            klass = info.get("CLASS") or _class_from_alleles(f[3], f[4], large_min)
            variants[(f[0], int(f[1]), f[3], f[4])] = klass
    return variants


def _fmt(value):
    return "NA" if isinstance(value, float) and math.isnan(value) else (
        f"{value:.4f}" if isinstance(value, float) else str(value)
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description="Sensitivity and precision of a call set against a truth set.")
    ap.add_argument("--truth", required=True)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--large-min", type=int, default=LARGE_INDEL_MIN)
    args = ap.parse_args(argv)

    stats = compare_variants(
        load_vcf(args.truth, args.large_min),
        load_vcf(args.calls, args.large_min),
        truth_svlen=load_svlens(args.truth),
        call_svlen=load_svlens(args.calls),
    )
    total = {k: sum(stats[c][k] for c in CLASSES) for k in ("truth", "tp", "fn", "calls", "fp")}
    total["sensitivity"] = _ratio(total["tp"], total["truth"])
    total["precision"] = _ratio(total["calls"] - total["fp"], total["calls"])
    columns = ("truth", "tp", "fn", "calls", "fp", "sensitivity", "precision")
    print("\t".join(("class",) + columns))
    for name, row in [*stats.items(), ("ALL", total)]:
        print("\t".join([name] + [_fmt(row[c]) for c in columns]))


if __name__ == "__main__":
    main()

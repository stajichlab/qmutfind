"""Compare a call set with a truth set: sensitivity and precision per variant class.

Both sets must be normalized first (`bcftools norm -f reference.fa`). A variant is
keyed by (contig, pos, ref, alt). A call matches only when all four fields are equal.
"""
import argparse
import math

from make_truth import _open_text

CLASSES = ("SNP", "SMALL_INDEL", "LARGE_INDEL")
LARGE_INDEL_MIN = 50


def _ratio(numerator, denominator):
    return numerator / denominator if denominator else float("nan")


def compare_variants(truth, calls):
    """Return {class: {truth, tp, fn, calls, fp, sensitivity, precision}}.
    `truth` and `calls` map a variant key to its class."""
    stats = {}
    for klass in CLASSES:
        truth_keys = {k for k, c in truth.items() if c == klass}
        call_keys = {k for k, c in calls.items() if c == klass}
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

    stats = compare_variants(load_vcf(args.truth, args.large_min), load_vcf(args.calls, args.large_min))
    total = {k: sum(stats[c][k] for c in CLASSES) for k in ("truth", "tp", "fn", "calls", "fp")}
    total["sensitivity"] = _ratio(total["tp"], total["truth"])
    total["precision"] = _ratio(total["calls"] - total["fp"], total["calls"])
    columns = ("truth", "tp", "fn", "calls", "fp", "sensitivity", "precision")
    print("\t".join(("class",) + columns))
    for name, row in [*stats.items(), ("ALL", total)]:
        print("\t".join([name] + [_fmt(row[c]) for c in columns]))


if __name__ == "__main__":
    main()

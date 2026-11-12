"""Score one run against its truth set and append the per-class counts to a TSV.

    python3 score_heldout_run.py --truth truth.norm.vcf --calls calls.norm.vcf --out results.tsv \
        --seed 3 --mode standard --tech art --coverage 30 [--contigs NC_1,NC_2]

Both VCFs must be normalized with `bcftools norm -f reference.fa`. Only truth variants and calls on
`--contigs` are scored, because the reads of a held-out run may come from a subset of the genome.
"""
import argparse
import os

from compare_calls import CLASSES, compare_variants, load_svlens, load_vcf

COLUMNS = ["seed", "mode", "tech", "coverage", "klass", "truth", "tp", "calls", "fp"]


def score_rows(truth_path, calls_path, contigs=None):
    truth = load_vcf(truth_path)
    calls = load_vcf(calls_path)
    call_svlen, truth_svlen = load_svlens(calls_path), load_svlens(truth_path)
    if contigs is not None:
        truth = {k: v for k, v in truth.items() if k[0] in contigs}
        calls = {k: v for k, v in calls.items() if k[0] in contigs}
    stats = compare_variants(truth, calls, truth_svlen=truth_svlen, call_svlen=call_svlen)
    return [
        {"klass": klass, "truth": stats[klass]["truth"], "tp": stats[klass]["tp"],
         "calls": stats[klass]["calls"], "fp": stats[klass]["fp"]}
        for klass in CLASSES
    ]


def main(argv=None):
    ap = argparse.ArgumentParser(description="Score one held-out run and append the rows to a TSV.")
    ap.add_argument("--truth", required=True)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--mode", required=True)
    ap.add_argument("--tech", required=True)
    ap.add_argument("--coverage", required=True)
    ap.add_argument("--contigs", help="comma-separated contig names to score (default: all)")
    args = ap.parse_args(argv)
    contigs = set(args.contigs.split(",")) if args.contigs else None
    new_file = not os.path.exists(args.out) or os.path.getsize(args.out) == 0
    with open(args.out, "a") as fh:
        if new_file:
            fh.write("\t".join(COLUMNS) + "\n")
        for row in score_rows(args.truth, args.calls, contigs):
            values = [args.seed, args.mode, args.tech, args.coverage, row["klass"], row["truth"], row["tp"], row["calls"], row["fp"]]
            fh.write("\t".join(str(v) for v in values) + "\n")


if __name__ == "__main__":
    main()

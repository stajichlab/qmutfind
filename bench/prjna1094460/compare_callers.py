"""Compare the call sets of five pipelines on the PRJNA1094460 samples (real reads, reference A1163). There is NO truth set.

The pipelines (per sample):
  rammap   qmutfind built-in caller on its own rammap alignments (what the browser tool does; reads mapped as single reads)
  bwa      qmutfind built-in caller on bwa-mem2 alignments (mates mapped as pairs)
  freebayes  qmutfind FreeBayes path on the same bwa-mem2 mapping
  bcf_raw  bwa-mem2 | bcftools mpileup | bcftools call --ploidy 1 -mv, all records
  bcf_flt  the same, at the qmutfind thresholds (depth 5, SNP fraction 0.8, indel fraction 0.5)

Every set is normalized in the same way (`bcftools norm -f ref -m -any`: multi-allelic records split, indels left-aligned) and a
call is (contig, position, REF, ALT). What the script measures is AGREEMENT and NOISE, not accuracy:
  - calls per pipeline, by class;
  - pairwise agreement (shared / only in one / Jaccard), for SNPs and for indels;
  - how many calls are supported by 1, 2, 3, 4 or 5 pipelines in the same sample;
  - private calls: a call that the same pipeline makes in no other sample. The samples share a CEA10 background (hundreds of
    differences from A1163), so most true calls are shared, and each sample has only a few tens of real new mutations;
  - a recall proxy: among (sample, variant) pairs where the variant is called by at least 3 pipelines in at least 10 of the
    13 samples ("background" variants, almost certainly real), the share that each pipeline calls. It depends on the pipelines
    that define the set, so it favours them; it is not sensitivity against a truth set.

Usage: python3 compare_callers.py --data DIR --ref plain.fa --out PREFIX [--samples a,b,c]
DIR holds calls/<name>.vcf (rammap) and compare/<name>.{bwa_builtin,bwa_freebayes}.vcf, .bcftools.{raw,filtered}.vcf
"""
import argparse
import collections
import itertools
import os
import subprocess
import sys

PIPELINES = ["rammap", "bwa", "freebayes", "bcf_raw", "bcf_flt"]


def path_of(data, name, pipeline):
    p = {
        "rammap": f"calls/{name}.vcf",
        "bwa": f"compare/{name}.bwa_builtin.vcf",
        "freebayes": f"compare/{name}.bwa_freebayes.vcf",
        "bcf_raw": f"compare/{name}.bcftools.raw.vcf",
        "bcf_flt": f"compare/{name}.bcftools.filtered.vcf",
    }[pipeline]
    return os.path.join(data, p)


def normalize(path, ref, cache_dir):
    """The path of the normalized copy (made once)."""
    out = os.path.join(cache_dir, os.path.basename(path) + ".norm.vcf")
    if not os.path.exists(out) or os.path.getmtime(out) < os.path.getmtime(path):
        subprocess.run(["bcftools", "norm", "-f", ref, "-m", "-any", path, "-o", out], check=True, stderr=subprocess.DEVNULL)
    return out


def klass(ref, alt):
    if alt.startswith("<"):
        return "SV"
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    return "INDEL"


def read_calls(path):
    """{(contig, pos, ref, alt): (af, dp)} from a normalized VCF. AF from INFO AF, or from FORMAT AD and DP (bcftools)."""
    calls = {}
    for line in open(path):
        if line[0] == "#":
            continue
        f = line.rstrip("\n").split("\t")
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        af = dp = None
        if "AF" in info and "DP" in info:
            af, dp = float(info["AF"].split(",")[0]), int(float(info["DP"]))
        elif len(f) > 9:
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            try:
                dp = int(fmt.get("DP") or info.get("DP"))
                ad = [int(x) for x in fmt["AD"].split(",")]
                af = ad[1] / dp if dp else None  # after `norm -m -any` each record has one ALT: AD is ref,alt
            except (KeyError, TypeError, ValueError, IndexError):
                pass
        calls[(f[0], int(f[1]), f[3], f[4])] = (af, dp)
    return calls


def jaccard(a, b):
    u = len(a | b)
    return len(a & b) / u if u else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--ref", required=True, help="the reference as plain FASTA with a .fai (for bcftools norm)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--samples", help="comma list (default: every sample in samples.tsv)")
    args = ap.parse_args()
    names = args.samples.split(",") if args.samples else [l.split("\t")[1] for l in open(os.path.join(args.data, "samples.tsv")).read().splitlines()[1:]]
    cache = os.path.join(args.data, "compare", "norm")
    os.makedirs(cache, exist_ok=True)

    calls = {}  # (sample, pipeline) -> {key: (af, dp)}
    usable = []
    for n in names:
        ok = True
        for p in PIPELINES:
            path = path_of(args.data, n, p)
            if not os.path.exists(path):
                ok = False
                break
        if not ok:
            print(f"skipping {n}: a call set is missing", file=sys.stderr)
            continue
        usable.append(n)
        for p in PIPELINES:
            calls[(n, p)] = read_calls(normalize(path_of(args.data, n, p), args.ref, cache))
    names = usable
    if len(names) < 2:
        sys.exit("need at least two samples with all five call sets")
    lines = []
    emit = lambda s="": (lines.append(s), print(s))

    emit(f"samples used: {len(names)}  ({', '.join(names)})")
    emit("\n== calls per pipeline (all samples pooled): SNP / INDEL / SV")
    for p in PIPELINES:
        c = collections.Counter(klass(k[2], k[3]) for n in names for k in calls[(n, p)])
        emit(f"  {p:10s} {c['SNP']:6d} {c['INDEL']:6d} {c['SV']:4d}   mean per sample {sum(c.values()) / len(names):.0f}")

    emit("\n== pairwise agreement within a sample (pooled over samples): shared, only first, only second, Jaccard")
    for cls in ("SNP", "INDEL"):
        emit(f"  -- {cls}")
        for a, b in itertools.combinations(PIPELINES, 2):
            both = oa = ob = 0
            for n in names:
                ka = {k for k in calls[(n, a)] if klass(k[2], k[3]) == cls}
                kb = {k for k in calls[(n, b)] if klass(k[2], k[3]) == cls}
                both += len(ka & kb); oa += len(ka - kb); ob += len(kb - ka)
            emit(f"  {a:10s} {b:10s} {both:6d} {oa:6d} {ob:6d}   {both / max(both + oa + ob, 1):.3f}")

    emit("\n== support: calls of the union, by the number of pipelines that make the call in the same sample")
    support = collections.Counter()
    by_pipeline_support = {p: collections.Counter() for p in PIPELINES}
    for n in names:
        union = set().union(*(calls[(n, p)].keys() for p in PIPELINES))
        for k in union:
            who = [p for p in PIPELINES if k in calls[(n, p)]]
            support[len(who)] += 1
            for p in who:
                by_pipeline_support[p][len(who)] += 1
    emit("  pipelines calling: " + "  ".join(f"{i}:{support[i]}" for i in range(1, 6)))
    emit("  per pipeline, how many of its calls are made by 1 / 2 / 3 / 4 / 5 pipelines:")
    for p in PIPELINES:
        emit(f"   {p:10s} " + "  ".join(f"{i}:{by_pipeline_support[p][i]:5d}" for i in range(1, 6)) + f"   made only by itself: {by_pipeline_support[p][1] / max(sum(by_pipeline_support[p].values()), 1):.1%}")

    emit("\n== private calls: made by this pipeline in one sample and in no other sample (pooled; per sample mean)")
    for p in PIPELINES:
        seen = collections.Counter(k for n in names for k in calls[(n, p)])
        priv = [sum(1 for k in calls[(n, p)] if seen[k] == 1) for n in names]
        tot = [len(calls[(n, p)]) for n in names]
        emit(f"  {p:10s} private {sum(priv):5d} of {sum(tot):6d} calls ({sum(priv) / max(sum(tot), 1):.1%}); per sample {sum(priv) / len(names):.1f}")

    emit("\n== recall proxy: background variants = called by >= 3 pipelines in >= 10 samples (or in all but 3 if fewer samples)")
    need = min(10, len(names) - 3) if len(names) >= 4 else len(names)
    count = collections.Counter()
    for n in names:
        for k in set().union(*(calls[(n, p)].keys() for p in PIPELINES)):
            if sum(1 for p in PIPELINES if k in calls[(n, p)]) >= 3:
                count[k] += 1
    background = {k for k, c in count.items() if c >= need}
    emit(f"  background variants: {len(background)} ({sum(1 for k in background if klass(k[2], k[3]) == 'SNP')} SNP, {sum(1 for k in background if klass(k[2], k[3]) == 'INDEL')} indel); (sample, variant) pairs: {len(background) * len(names)}")
    for p in PIPELINES:
        hit = {c: 0 for c in ("SNP", "INDEL")}
        tot = {c: 0 for c in ("SNP", "INDEL")}
        for n in names:
            for k in background:
                c = klass(k[2], k[3])
                if c in tot:
                    tot[c] += 1
                    hit[c] += k in calls[(n, p)]
        emit(f"  {p:10s} SNP {hit['SNP']:5d}/{tot['SNP']:5d} = {hit['SNP'] / max(tot['SNP'], 1):.3f}   INDEL {hit['INDEL']:4d}/{tot['INDEL']:4d} = {hit['INDEL'] / max(tot['INDEL'], 1):.3f}")

    emit("\n== calls made by only one pipeline (pooled): number, median AF, median depth, share with AF < 0.5")
    for p in PIPELINES:
        afs, dps = [], []
        for n in names:
            for k, (af, dp) in calls[(n, p)].items():
                if not any(k in calls[(n, q)] for q in PIPELINES if q != p):
                    if af is not None: afs.append(af)
                    if dp is not None: dps.append(dp)
        med = lambda v: sorted(v)[len(v) // 2] if v else float("nan")
        emit(f"  {p:10s} {len(afs):5d}   AF {med(afs):.2f}   DP {med(dps):5.1f}   AF<0.5: {sum(1 for a in afs if a < 0.5) / max(len(afs), 1):.1%}")

    with open(args.out + ".txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    # one row per sample and pipeline
    with open(args.out + ".per_sample.tsv", "w") as fh:
        fh.write("sample\tpipeline\tsnp\tindel\tsv\n")
        for n in names:
            for p in PIPELINES:
                c = collections.Counter(klass(k[2], k[3]) for k in calls[(n, p)])
                fh.write(f"{n}\t{p}\t{c['SNP']}\t{c['INDEL']}\t{c['SV']}\n")


if __name__ == "__main__":
    main()

"""Does pair-aware mapping change the calls? Compares rammap single-read calls (calls_single/) with rammap paired calls
(calls_paired/), both from the same binary, against the bwa-mem2 pipelines on the PRJNA1094460 samples. There is NO truth set.

Reported (pooled over samples):
  - calls per sample, SNP / INDEL;
  - calls gained and lost by pairing, and how many of those the other pipelines (bwa built-in, FreeBayes, bcftools raw) make in the
    same sample (agreement is support, not proof);
  - the recall proxy of compare_callers.py (background variants called by >= 3 of the five earlier pipelines in >= 10 samples);
  - private calls (made in one sample only);
  - Jaccard of each rammap set with the bwa built-in set.
Usage: python3 compare_pairing.py --data DIR --ref plain.fa --out PREFIX
"""
import argparse, collections, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_callers as cc

OTHERS = ["bwa", "freebayes", "bcf_raw", "bcf_flt"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    names = [l.split("\t")[1] for l in open(os.path.join(a.data, "samples.tsv")).read().splitlines()[1:]]
    cache = os.path.join(a.data, "compare", "norm")
    os.makedirs(cache, exist_ok=True)
    paths = {"single": "calls_single/{n}.vcf", "paired": "calls_paired/{n}.vcf"}
    calls = {}
    for n in names:
        for p in OTHERS:
            calls[(n, p)] = cc.read_calls(cc.normalize(cc.path_of(a.data, n, p), a.ref, cache))
        for p, t in paths.items():
            src = os.path.join(a.data, t.format(n=n))
            # a distinct cache name per set: the file names are the same in both directories
            dst = os.path.join(cache, f"{p}.{n}.vcf")
            import subprocess
            subprocess.run(["bcftools", "norm", "-f", a.ref, "-m", "-any", src, "-o", dst], check=True, stderr=subprocess.DEVNULL)
            calls[(n, p)] = cc.read_calls(dst)
    lines = []
    emit = lambda s="": (lines.append(s), print(s))
    k = lambda key: cc.klass(key[2], key[3])

    emit(f"samples: {len(names)}")
    emit("\n== calls (pooled over samples) SNP / INDEL, mean per sample")
    for p in ["single", "paired"] + OTHERS:
        c = collections.Counter(k(x) for n in names for x in calls[(n, p)])
        emit(f"  {p:10s} {c['SNP']:5d} {c['INDEL']:4d}   {sum(c.values()) / len(names):.0f}")

    emit("\n== calls gained / lost by pairing (pooled), and support from other pipelines in the same sample")
    for label, src, dst in (("gained", "single", "paired"), ("lost", "paired", "single")):
        tot = collections.Counter(); sup = collections.Counter(); anyp = collections.Counter(); bwa_only = collections.Counter()
        for n in names:
            for x in set(calls[(n, dst)]) - set(calls[(n, src)]):  # in dst and not in src
                c = k(x); tot[c] += 1
                who = [p for p in OTHERS if x in calls[(n, p)]]
                anyp[c] += bool(who)
                sup[c] += len(who) >= 2
        for c in ("SNP", "INDEL"):
            emit(f"  {label:7s} {c:5s} {tot[c]:4d}   made by >=1 other pipeline: {anyp[c]:4d} ({anyp[c] / max(tot[c], 1):.0%})   by >=2: {sup[c]:4d} ({sup[c] / max(tot[c], 1):.0%})")

    emit("\n== recall proxy (background variants as in compare_callers.py)")
    need = min(10, len(names) - 3)
    five = ["bwa", "freebayes", "bcf_raw", "bcf_flt"]  # the earlier single-read rammap set is added from calls_single
    count = collections.Counter()
    for n in names:
        sets = {p: calls[(n, p)] for p in five}
        sets["single"] = calls[(n, "single")]
        for x in set().union(*sets.values()):
            if sum(x in s for s in sets.values()) >= 3:
                count[x] += 1
    background = {x for x, c in count.items() if c >= need}
    emit(f"  background variants: {len(background)} (defined by single, bwa, freebayes, bcf_raw, bcf_flt; this favours 'single', not 'paired')")
    for p in ["single", "paired"] + OTHERS:
        hit = collections.Counter(); tot = collections.Counter()
        for n in names:
            for x in background:
                tot[k(x)] += 1; hit[k(x)] += x in calls[(n, p)]
        emit(f"  {p:10s} SNP {hit['SNP']}/{tot['SNP']} = {hit['SNP'] / max(tot['SNP'], 1):.3f}   INDEL {hit['INDEL']}/{tot['INDEL']} = {hit['INDEL'] / max(tot['INDEL'], 1):.3f}")

    emit("\n== private calls (this set, one sample only)")
    for p in ["single", "paired", "bwa"]:
        seen = collections.Counter(x for n in names for x in calls[(n, p)])
        priv = sum(1 for n in names for x in calls[(n, p)] if seen[x] == 1)
        tot = sum(len(calls[(n, p)]) for n in names)
        emit(f"  {p:10s} {priv:5d} of {tot:6d} ({priv / tot:.1%}); per sample {priv / len(names):.1f}")

    emit("\n== Jaccard with the bwa built-in set (pooled) and calls only in bwa")
    for p in ["single", "paired"]:
        for cls in ("SNP", "INDEL"):
            both = oa = ob = 0
            for n in names:
                x = {y for y in calls[(n, p)] if k(y) == cls}; b = {y for y in calls[(n, "bwa")] if k(y) == cls}
                both += len(x & b); oa += len(x - b); ob += len(b - x)
            emit(f"  {p:7s} {cls:5s} shared {both:5d}  only {p} {oa:4d}  only bwa {ob:4d}  Jaccard {both / (both + oa + ob):.3f}")
    open(a.out + ".txt", "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

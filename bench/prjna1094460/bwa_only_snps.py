"""Why does bwa-mem2 + built-in caller call SNPs that paired rammap does not? For each SNP in the bwa call set and not in the
rammap paired call set of the same sample, look at what rammap saw at that site (its evidence VCF: sites with >= 2 alt reads
at relaxed thresholds, and the rammap depth in the depth BED). Pooled over samples. Usage: bwa_only_snps.py DATA OUT.tsv"""
import collections, os, sys, bisect
data, out = sys.argv[1], sys.argv[2]
names = [l.split("\t")[1] for l in open(os.path.join(data, "samples.tsv")).read().splitlines()[1:]]

def snps(path):
    d = {}
    for l in open(path):
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t")
        if len(f[3]) == 1 and len(f[4]) == 1:
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            d[(f[0], int(f[1]), f[3], f[4])] = (float(info["AF"]), int(info["DP"]), int(info["AO"]))
    return d

def bed(path):
    d = collections.defaultdict(list)
    for l in open(path):
        c, s, e, m = l.split(); d[c].append((int(s), int(e), float(m)))
    return d

def depth_at(b, c, pos):
    for s, e, m in b[c]:
        if s < pos <= e: return m
    return None

rows = []
cats = collections.Counter()
for n in names:
    bwa = snps(f"{data}/compare/norm/{n}.bwa_builtin.vcf.norm.vcf")
    ram = snps(f"{data}/compare/norm/paired.{n}.vcf")
    evi = snps(f"{data}/calls_paired/{n}.evidence.vcf")
    bd = bed(f"{data}/calls_paired/{n}.depth.bed")
    recur = collections.Counter()
    for m in names:
        for k in snps(f"{data}/compare/norm/{m}.bwa_builtin.vcf.norm.vcf"): recur[k] += 1
    for k, (af, dp, ao) in bwa.items():
        if k in ram: continue
        e = evi.get(k)
        if e is None:
            # same position, other alt? or no evidence at all
            other = [x for x in evi if x[:3] == k[:3]]
            cat = "rammap evidence: other ALT at site" if other else "rammap: no alt evidence (<2 alt reads at mapq>=20)"
        else:
            eaf, edp, eao = e
            if eaf >= 0.8 and edp < 5: cat = "rammap saw it, depth < 5"
            elif eaf >= 0.8: cat = "rammap saw it at AF>=0.8 (filtered otherwise: QUAL/depth ratio)"
            else: cat = "rammap saw it, AF < 0.8"
        cats[cat] += 1
        rows.append((n, *k, af, dp, ao, cat, e[0] if e else "", e[1] if e else "", e[2] if e else "", depth_at(bd, k[0], k[1]) or "", recur[k]))
with open(out, "w") as fh:
    fh.write("sample\tcontig\tpos\tref\talt\tbwa_af\tbwa_dp\tbwa_ao\tcategory\trammap_af\trammap_dp\trammap_ao\trammap_window_depth\tbwa_samples\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
tot = sum(cats.values())
print("bwa-only SNPs:", tot)
for c, v in cats.most_common(): print(f"  {v:5d} {v/tot:5.1%}  {c}")

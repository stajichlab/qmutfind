"""Count introns derived as gaps between consecutive exons of one mRNA, and check GT..AG ends.
Usage: python3 -I derive_intron_count.py <annotation.gff.gz> <genome.fa.gz>"""
import gzip, re, sys
from collections import defaultdict
gff, fa = sys.argv[1], sys.argv[2]
op = lambda p: gzip.open(p, "rt") if p.endswith(".gz") else open(p)
seq = {}; name = None; buf = []
with op(fa) as fh:
    for l in fh:
        if l.startswith(">"):
            if name: seq[name] = "".join(buf).upper()
            name = l[1:].split()[0]; buf = []
        else: buf.append(l.strip())
    if name: seq[name] = "".join(buf).upper()
comp = str.maketrans("ACGT", "TGCA")
mrna = {}; exons = defaultdict(list)
def attr(s, k):
    m = re.search(r"(?:^|;)" + k + r"=([^;]+)", s); return m.group(1) if m else None
with op(gff) as fh:
    for l in fh:
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9: continue
        if f[2] == "mRNA": mrna[attr(f[8], "ID")] = (f[0], f[6])
        elif f[2] == "exon":
            for p in (attr(f[8], "Parent") or "").split(","): exons[p].append((int(f[3]), int(f[4])))
n = 0; lens = []; ends = defaultdict(int); multi = 0; overlap = 0
for mid, (chrom, strand) in mrna.items():
    ex = sorted(exons.get(mid, []))
    if len(ex) > 1: multi += 1
    for (a1, b1), (a2, b2) in zip(ex, ex[1:]):
        s, e = b1 + 1, a2 - 1
        if e < s: overlap += 1; continue
        n += 1; lens.append(e - s + 1)
        if chrom in seq:
            ig = seq[chrom][s - 1:e]
            if strand == "-": ig = ig.translate(comp)[::-1]
            ends[ig[:2] + ".." + ig[-2:]] += 1
lens.sort()
print(f"mRNA: {len(mrna)}; multi-exon mRNA: {multi}")
print(f"derived introns: {n}; adjacent exons with no gap: {overlap}")
print(f"intron length min/median/max: {lens[0]}/{lens[len(lens)//2]}/{lens[-1]}")
tot = sum(ends.values())
for k, v in sorted(ends.items(), key=lambda x: -x[1])[:6]: print(f"  {k}: {v} ({100*v/tot:.1f}%)")

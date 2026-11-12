"""Count mRNAs with a derivable 5' or 3' UTR (exon bases outside the CDS span).
Usage: python3 -I derive_utr_count.py <annotation.gff.gz>
"""
import gzip
import re
import sys
from collections import defaultdict

path = sys.argv[1]
opener = gzip.open if path.endswith(".gz") else open
mrna = {}
exons = defaultdict(list)
cds = defaultdict(list)


def attr(s, key):
    m = re.search(r"(?:^|;)" + key + r"=([^;]+)", s)
    return m.group(1) if m else None


with opener(path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        typ, start, end, strand, a = f[2], int(f[3]), int(f[4]), f[6], f[8]
        if typ == "mRNA":
            mrna[attr(a, "ID")] = (f[0], start, end, strand)
        elif typ == "exon":
            for p in (attr(a, "Parent") or "").split(","):
                exons[p].append((start, end))
        elif typ == "CDS":
            for p in (attr(a, "Parent") or "").split(","):
                cds[p].append((start, end))

total = len(mrna)
no_cds = no_exon = 0
five = three = both = none = 0
five_bp = []
three_bp = []
for mid, (chrom, s, e, strand) in mrna.items():
    if not cds.get(mid):
        no_cds += 1
        continue
    if not exons.get(mid):
        no_exon += 1
        continue
    cmin = min(a for a, _ in cds[mid])
    cmax = max(b for _, b in cds[mid])
    left = sum(max(0, min(b, cmin - 1) - a + 1) for a, b in exons[mid])
    right = sum(max(0, b - max(a, cmax + 1) + 1) for a, b in exons[mid])
    f5, f3 = (left, right) if strand == "+" else (right, left)
    if f5 > 0:
        five_bp.append(f5)
    if f3 > 0:
        three_bp.append(f3)
    if f5 > 0 and f3 > 0:
        both += 1
    elif f5 > 0:
        five += 1
    elif f3 > 0:
        three += 1
    else:
        none += 1


def med(x):
    x = sorted(x)
    return x[len(x) // 2] if x else 0


print(f"mRNA total: {total}")
print(f"mRNA with no CDS: {no_cds}; with no exon: {no_exon}")
print(f"both UTRs: {both}; only 5': {five}; only 3': {three}; neither: {none}")
print(f"5' UTR present: {len(five_bp)} (median {med(five_bp)} bp); 3' UTR present: {len(three_bp)} (median {med(three_bp)} bp)")

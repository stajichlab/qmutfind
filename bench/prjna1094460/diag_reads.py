"""For the SNPs that bwa-mem2 + built-in caller call in WT30 and paired rammap does not (bwa_only_snps.tsv), take the reads that
carry the ALT base in the bwa-mem2 BAM and say what the rammap single-read alignment of the same read looks like.
Usage: diag_reads.py DIAG_DIR BWA_ONLY.tsv OUT.tsv"""
import collections, gzip, re, subprocess, sys
diag, sites_tsv, out = sys.argv[1:4]
sites = [l.rstrip("\n").split("\t") for l in open(sites_tsv)][1:]
sites = [s for s in sites if s[0] == "WT30"]
need = {}  # (readname, mate) -> None, filled later
site_reads = []  # (site index, readname, mate, flag, mapq, cigar, qpos, nm)
cig = re.compile(r"(\d+)([MIDNSHP=X])")

def base_at(pos1, start1, cigar, seq):
    """(base, query index, clipped_left, clipped_right, indel_in_read) at 1-based ref pos, or None."""
    r, q = start1, 0
    ops = cig.findall(cigar)
    for i, (n, o) in enumerate(ops):
        n = int(n)
        if o in "M=X":
            if r <= pos1 < r + n:
                return seq[q + pos1 - r], q + pos1 - r
            r += n; q += n
        elif o in "IS": q += n
        elif o in "DN":
            if r <= pos1 < r + n: return None
            r += n
    return None

for i, s in enumerate(sites):
    c, pos, ref, alt = s[1], int(s[2]), s[3], s[4]
    res = subprocess.run(["samtools", "view", f"{diag}/WT30.bwa.bam", f"{c}:{pos}-{pos}"], capture_output=True, text=True, check=True).stdout
    for l in res.splitlines():
        f = l.split("\t")
        flag = int(f[1])
        if flag & 0x904: continue
        b = base_at(pos, int(f[3]), f[5], f[9])
        if b is None or b[0] != alt: continue
        nm = next((int(x[5:]) for x in f[11:] if x.startswith("NM:i:")), -1)
        mate = 1 if flag & 0x40 else 2
        site_reads.append((i, f[0], mate, flag, int(f[4]), f[5], b[1], nm, len(f[9])))
        need[(f[0], mate)] = None

for m in (1, 2):
    for l in gzip.open(f"{diag}/WT30.rammap.m{m}.tsv.gz", "rt"):
        f = l.rstrip("\n").split("\t")
        k = (f[0], m)
        if k in need: need[k] = f

cats = collections.Counter(); percat_sites = collections.defaultdict(collections.Counter)
rows = []
for (i, name, mate, flag, mq, cigar, qi, nm, rl) in site_reads:
    s = sites[i]; c, pos = s[1], int(s[2])
    r = need[(name, mate)]
    if r is None: cat = "read not in rammap output"
    elif r[3] == "*": cat = "rammap: unmapped"
    else:
        rc, rpos, rcig, rmq = r[3], int(r[4]), r[6], int(r[7])
        if rc != c or abs(rpos - int(pos)) > 2000: cat = "rammap: mapped elsewhere"
        elif rmq < 20: cat = "rammap: same place, mapq < 20"
        else:
            # same place, mapq ok: clipped or alt base lost?
            clipped = "S" in cigar
            rc_ops = cig.findall(rcig); covered_len = sum(int(n) for n, o in rc_ops if o in "MIS")  # rammap emits M/I/D only; length of query aligned
            if sum(int(n) for n, o in rc_ops if o in "MI") < rl - 0: cat = "rammap: same place, read clipped (aligned < read length)"
            elif "I" in rcig or "D" in rcig: cat = "rammap: same place, alignment has an indel"
            else: cat = "rammap: same place, full-length gap-free alignment (alt base should count)"
    cats[cat] += 1
    percat_sites[i][cat] += 1
    rows.append((s[0], c, pos, i, name, mate, mq, cigar, qi, nm, cat, "" if r is None else "\t".join(r[3:8])))
with open(out, "w") as fh:
    fh.write("sample\tcontig\tpos\tsite\tread\tmate\tbwa_mapq\tbwa_cigar\tbwa_qpos\tbwa_NM\tcategory\trammap_contig\trammap_pos\trammap_strand\trammap_cigar\trammap_mapq\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
tot = sum(cats.values())
print(f"sites {len(sites)}; alt-carrying bwa reads {tot}")
for c, v in cats.most_common(): print(f"  {v:6d} {v/tot:5.1%} {c}")
# the dominant reason per site (the category of most of its reads)
dom = collections.Counter(max(v, key=v.get) for v in percat_sites.values())
print("dominant category per site:", dict(dom), " sites with no alt read found in BAM:", len(sites) - len(percat_sites))

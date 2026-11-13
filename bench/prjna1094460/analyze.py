import sys, statistics, collections, urllib.parse, gzip
P = "/bigdata/stajichlab/jstajich/projects/qmutfind_data/prjna1094460"
C = P + "/calls"
names = ["WT30","WT_A29","WT_B23","WT_C30","WT_D30","crpA_control","crpA_B23","crpA_C10","crpA_D15","aceA_control","aceA_B22","aceA_C22","aceA_D22"]
group = lambda n: "WT" if n.startswith("WT") else ("crpA" if n.startswith("crpA") else "aceA")

# --- 1. deletion strains: windows with no coverage in every sample of one genotype but covered in the others
depth = {}
for n in names:
    d = {}
    for l in open(f"{C}/{n}.depth.bed"):
        c, s, e, v = l.split("\t"); d[(c, int(s))] = (int(e), float(v))
    depth[n] = d
keys = list(depth["WT30"])
for g in ("crpA", "aceA"):
    mine = [n for n in names if group(n) == g]; others = [n for n in names if group(n) != g]
    hits = [k for k in keys if all(depth[n][k][1] < 2 for n in mine) and statistics.median(depth[n][k][1] for n in others) >= 10]
    print(f"\n[{g}] windows with depth < 2 in all {len(mine)} {g} samples and median >= 10 in the other {len(others)}: {len(hits)}")
    # merge adjacent windows
    runs = []
    for k in sorted(hits):
        if runs and runs[-1][0] == k[0] and runs[-1][2] == k[1]: runs[-1][2] = depth["WT30"][k][0]
        else: runs.append([k[0], k[1], depth["WT30"][k][0]])
    for c, s, e in runs[:10]: print("   ", c, s, e, f"({e-s} bp)")

# crpA candidate gene AFUB_069550 at DS499598.1:1871265-1874828
def window_mean(n, c, s, e):
    vals = [v for (cc, ss), (ee, v) in depth[n].items() if cc == c and ss < e and ee > s]
    return statistics.mean(vals) if vals else float("nan")
print("\nmean depth over the copper-transporting ATPase candidate AFUB_069550 (DS499598.1:1871265-1874828):")
for n in names: print(f"   {n:14s} {window_mean(n,'DS499598.1',1871265,1874828):7.1f}   genome median {statistics.median(v for _, v in depth[n].values()):.1f}")

# --- 2. annotation
ann = {}
for n in names:
    a = {}
    for l in open(f"{C}/{n}.annot.vcf"):
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t"); info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        a[(f[0], int(f[1]), f[3], f[4].split(",")[0])] = info
    ann[n] = a

# --- 3. the combined table
rows = [l.rstrip("\n").split("\t") for l in open(C + "/all.combined.tsv")]
hdr = rows[0]; k = len(names); first = 6
status_col = {n: hdr.index(f"{n}.status") for n in names}
cell_col = {n: hdr.index(n) for n in names}
new_by_sample = collections.defaultdict(list)
for r in rows[1:]:
    key = (r[0], int(r[1]), r[2], r[3])
    for n in names[1:]:
        if r[status_col[n]] == "new": new_by_sample[n].append(key)

def info_for(key, samples):
    for n in samples:
        if key in ann[n]: return ann[n][key]
    return {}
print("\n[new vs WT30] counts and effects per sample")
for n in names[1:]:
    eff = collections.Counter(info_for(k_, [n]).get("EFFECT", "?").split(",")[0] for k_ in new_by_sample[n])
    print(f"   {n:14s} {len(new_by_sample[n]):3d}  " + ", ".join(f"{e}:{c}" for e, c in eff.most_common(6)))

# recurrent: the same variant new in >= 2 samples, and the same gene hit in >= 2 samples
by_variant = collections.defaultdict(list)
for n, lst in new_by_sample.items():
    for key in lst: by_variant[key].append(n)
print("\n[same variant new in 2+ samples]")
for key, ns in sorted(by_variant.items(), key=lambda x: -len(x[1])):
    if len(ns) >= 2:
        i = info_for(key, ns); print(f"   {key[0]}:{key[1]} {key[2]}>{key[3]}  {len(ns)} samples {','.join(ns)}  {i.get('EFFECT','?')} {i.get('GENE','')} {urllib.parse.unquote(i.get('PRODUCT',''))[:50]}")
gene_hits = collections.defaultdict(set)
for key, ns in by_variant.items():
    i = info_for(key, ns)
    if i.get("EFFECT", "").split(",")[0] in ("missense_variant","frameshift_variant","stop_gained","start_lost","stop_lost","inframe_deletion","inframe_insertion","splice_donor_variant","splice_acceptor_variant", "nonsense", "missense", "frameshift"):
        gene_hits[i.get("GENE", "?")].update(ns)
print("\n[genes with a coding change new in 2+ different samples]")
for g, ns in sorted(gene_hits.items(), key=lambda x: -len(x[1])):
    if len(ns) >= 2 and g not in ("", "?", "."): print("  ", g, len(ns), ",".join(sorted(ns)))
print("\neffect labels seen:", sorted({i.get("EFFECT","?").split(",")[0] for a in ann.values() for i in a.values()}))

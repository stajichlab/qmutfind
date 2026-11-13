"""A deterministic test set of variants for comparing annotators: python3 make_variants.py GFF REF.fa SEED OUT.vcf
  - every possible SNV at every position of 150 random protein-coding genes and 300 bp on each side;
  - every possible SNV within 12 bp of each side of every exon boundary of 600 other random genes (splice sites, UTR and intron edges);
  - 4,000 random insertions and deletions of 1-9 bp inside genes (and 300 bp around them), and 3,000 near exon boundaries.
Only the eight nuclear chromosomes (NC_007194.1 to NC_007201.1). The VCF has no genotype columns. Normalize it with `bcftools norm -f REF` before use."""
import random, sys, collections
gff, ref, seed, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
rng = random.Random(seed)
NUC = {f"NC_0071{n}.1" for n in range(94, 102)}
genome = {}
name = None
for line in open(ref):
    if line[0] == ">":
        name = line[1:].split()[0]; genome[name] = []
    else:
        genome[name].append(line.strip().upper())
genome = {k: "".join(v) for k, v in genome.items() if k in NUC}

genes, mrna_gene, exons = {}, {}, collections.defaultdict(list)
for line in open(gff):
    if line[0] == "#": continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[0] not in NUC: continue
    attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    if f[2] == "gene" and attrs.get("gene_biotype") == "protein_coding":
        genes[attrs["ID"]] = (f[0], int(f[3]), int(f[4]))
    elif f[2] == "mRNA":
        mrna_gene[attrs["ID"]] = attrs.get("Parent")
    elif f[2] == "exon":
        for p in attrs.get("Parent", "").split(","):
            exons[p].append((f[0], int(f[3]), int(f[4])))
gene_exons = collections.defaultdict(set)
for m, lst in exons.items():
    g = mrna_gene.get(m)
    if g in genes:
        gene_exons[g].update(lst)
ids = sorted(g for g in genes if g in gene_exons)
rng.shuffle(ids)
full, edge = ids[:150], ids[150:750]

variants = set()  # (chrom, pos, ref, alt)
def snvs(chrom, lo, hi):
    seq = genome[chrom]
    for pos in range(max(1, lo), min(len(seq), hi) + 1):
        r = seq[pos - 1]
        if r not in "ACGT": continue
        for a in "ACGT":
            if a != r: variants.add((chrom, pos, r, a))
for g in full:
    c, s, e = genes[g]; snvs(c, s - 300, e + 300)
for g in edge:
    for c, s, e in gene_exons[g]:
        snvs(c, s - 12, s + 12); snvs(c, e - 12, e + 12)

def indel(chrom, pos):
    seq = genome[chrom]
    if pos < 2 or pos > len(seq) - 20: return
    n = rng.randint(1, 9)
    if rng.random() < 0.5:  # deletion of n bases after the anchor at pos
        r = seq[pos - 1: pos + n]
        if set(r) <= set("ACGT"): variants.add((chrom, pos, r, r[0]))
    else:  # insertion of n random bases after the anchor
        anchor = seq[pos - 1]
        if anchor in "ACGT": variants.add((chrom, pos, anchor, anchor + "".join(rng.choice("ACGT") for _ in range(n))))
for _ in range(4000):
    g = rng.choice(full); c, s, e = genes[g]; indel(c, rng.randint(s - 300, e + 300))
edges = [(c, s, e) for g in full + edge for c, s, e in gene_exons[g]]
for _ in range(3000):
    c, s, e = rng.choice(edges); indel(c, (s if rng.random() < 0.5 else e) + rng.randint(-12, 12))

order = {c: i for i, c in enumerate(sorted(genome))}
with open(out, "w") as fh:
    fh.write("##fileformat=VCFv4.2\n")
    for c in sorted(genome): fh.write(f"##contig=<ID={c},length={len(genome[c])}>\n")
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for c, p, r, a in sorted(variants, key=lambda v: (order[v[0]], v[1], v[2], v[3])):
        fh.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\n")
print(f"{len(variants)} variants in {len(full)} full genes and {len(edge)} edge genes", file=sys.stderr)

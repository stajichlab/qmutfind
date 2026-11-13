"""Compare `qmutfind annotate` with snpEff on the same VCF: python3 compare_snpeff.py QMUTFIND.vcf SNPEFF.vcf OUT_PREFIX
snpEff must be run with -canon (longest-CDS transcript, as qmutfind uses) and -ud 0 (no upstream/downstream labels, so a site outside genes is
intergenic_region). This is NOT a test against a truth set: both tools can be wrong, and they define some classes differently.
Every variant ends in one of three outcomes: AGREE (snpEff gives a term that means the qmutfind label), CONVENTION (they differ in a known definition,
see `convention`), DIFFER (anything else). Writes OUT_PREFIX.report.txt, OUT_PREFIX.disagreements.tsv (CONVENTION and DIFFER rows) and
OUT_PREFIX.aa_disagreements.tsv."""
import collections, re, sys

# qmutfind EFFECT label -> the snpEff (Sequence Ontology) terms that mean the same. A snpEff annotation can hold several terms joined by '&'.
ACCEPT = {
    "synonymous": {"synonymous_variant", "stop_retained_variant", "start_retained_variant"},
    "missense": {"missense_variant"},
    "nonsense": {"stop_gained"},
    "start_lost": {"start_lost"},
    "stop_lost": {"stop_lost"},
    "frameshift": {"frameshift_variant"},
    "inframe_insertion": {"conservative_inframe_insertion", "disruptive_inframe_insertion", "inframe_insertion"},
    "inframe_deletion": {"conservative_inframe_deletion", "disruptive_inframe_deletion", "inframe_deletion"},
    "splice_site": {"splice_acceptor_variant", "splice_donor_variant"},
    "5UTR_derived": {"5_prime_UTR_variant", "5_prime_UTR_premature_start_codon_gain_variant"},
    "3UTR_derived": {"3_prime_UTR_variant"},
    "intron_derived": {"intron_variant"},
    "non_coding": {"non_coding_transcript_exon_variant", "non_coding_transcript_variant"},
    "intergenic": {"intergenic_region"},
    "large_indel_in_cds": {"frameshift_variant", "conservative_inframe_deletion", "disruptive_inframe_deletion", "conservative_inframe_insertion",
                           "disruptive_inframe_insertion", "stop_gained", "start_lost", "stop_lost", "exon_loss_variant", "feature_ablation", "transcript_ablation"},
    "complex": set(),  # a variant over an exon boundary: judged by hand
}

def convention(label, terms):
    """The known difference of definition that explains a disagreement, or None."""
    if label == "start_lost" and "initiator_codon_variant" in terms:
        return "qmutfind start_lost / snpEff initiator_codon_variant (a change to another start codon of the table, e.g. ATG to CTG)"
    if terms == {"splice_region_variant"}:
        return "snpEff gives only splice_region_variant (its convention near an exon edge); qmutfind labels the position by the exons"
    if label == "non_coding" and terms & {"intragenic_variant", "intron_variant", "splice_region_variant", "non_coding_transcript_exon_variant"}:
        return "qmutfind non_coding / snpEff intragenic_variant (a feature that is not a protein-coding transcript)"
    return None

THREE = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
         "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U", "Pyl": "O"}

def one_letter(hgvs_p):
    """p.Val29Ile -> V29I ; p.Glu3Aspfs*4 -> E3Dfs*4 ; p.Met1? -> M1? ; '' when it cannot be read."""
    m = re.fullmatch(r"p\.\(?([A-Z][a-z]{2})(\d+)(.*?)\)?", hgvs_p)
    if not m: return ""
    ref, pos, rest = m.groups()
    rest = re.sub(r"[A-Z][a-z]{2}", lambda x: THREE.get(x.group(0), x.group(0)), rest)
    return f"{THREE.get(ref, ref)}{pos}{rest}"

def ref_and_position(aa):
    m = re.match(r"([A-Z*])(\d+)", aa)
    return (m.group(1), m.group(2)) if m else None

def info_of(field):
    d = {}
    for kv in field.split(";"):
        k, _, v = kv.partition("=")
        d[k] = v
    return d

def read(path):
    out = {}
    for line in open(path):
        if line[0] == "#": continue
        f = line.rstrip("\n").split("\t")
        out[(f[0], int(f[1]), f[3].upper(), f[4].upper())] = info_of(f[7])  # snpEff upper-cases the bases of a soft-masked reference
    return out

klass = lambda ref, alt: "SNV" if len(ref) == 1 and len(alt) == 1 else "indel"

def main():
    q, s = read(sys.argv[1]), read(sys.argv[2])
    assert q.keys() == s.keys(), "the two VCFs hold different variants"
    total = collections.Counter(); agree = collections.Counter(); conv = collections.Counter(); differ = collections.Counter()
    why = collections.Counter(); first_term = collections.Counter()
    rows = []
    aa_total = collections.Counter(); aa_equal = collections.Counter(); aa_note = collections.Counter(); aa_bad = []
    gene_total = collections.Counter(); gene_same = collections.Counter(); gene_explained = collections.Counter(); gene_bad = []
    for key, qi in q.items():
        ref, alt = key[2], key[3]
        c = klass(ref, alt)
        label = qi.get("EFFECT", "").split(",")[0]
        entries = [e.split("|") for e in s[key].get("ANN", "").split(",") if e]
        if not entries: continue
        qtx = qi.get("TRANSCRIPT", "").removeprefix("rna-")
        entry = next((e for e in entries if len(e) > 6 and e[6] == qtx), entries[0])  # the snpEff entry of the same transcript, else its first
        terms = set(entry[1].split("&"))
        total[(c, label)] += 1
        ok = bool(terms & ACCEPT.get(label, set()))
        if ok:
            agree[(c, label)] += 1
        else:
            reason = convention(label, terms)
            if reason:
                conv[(c, label)] += 1; why[reason] += 1
            else:
                differ[(c, label)] += 1; first_term[(label, entry[1])] += 1
            rows.append((*key, label, entry[1], qi.get("GENE", ""), entry[4], qi.get("TRANSCRIPT", ""), entry[6], qi.get("AA", ""), entry[10], "convention" if reason else "differ"))
        # gene and transcript identity
        if label != "intergenic" and qtx:
            gene_total[c] += 1
            if entry[6] == qtx:
                gene_same[c] += 1
            elif any(e[4].removeprefix("gene-") == qi.get("GENE", "") for e in entries):
                gene_explained[c] += 1  # snpEff reports the gene too, but its canonical transcript or its first entry differs (overlapping genes, equal CDS lengths)
            elif len(gene_bad) < 200: gene_bad.append((*key, label, qi.get("GENE", ""), qtx, entry[4], entry[6], entry[1]))
        # protein change, where the classes agree and both give one
        qa = qi.get("AA", "").split(",")[0]
        sa = one_letter(entry[10]) if len(entry) > 10 else ""
        if ok and qa and qa != "." and sa and not sa.endswith("="):
            aa_total[(c, label)] += 1
            if label in ("start_lost", "stop_lost", "frameshift", "large_indel_in_cds"):
                # the tools write these differently (qmutfind M1V / *244R / E3Dfs*4; snpEff p.Met1? / p.Ter244Argext*? / p.Glu3fs): the reference residue and position are compared
                same = ref_and_position(qa) == ref_and_position(sa)
                aa_note[(c, label)] += 1
            else:
                same = qa == sa
            if same: aa_equal[(c, label)] += 1
            elif len(aa_bad) < 4000: aa_bad.append((*key, label, qa, sa, entry[10]))

    L = []
    emit = lambda t="": (L.append(t), print(t))
    num = lambda d, c: sum(v for (cc, _), v in d.items() if cc == c)
    emit(f"variants compared: {sum(total.values()):,}   (SNV {num(total, 'SNV'):,}, indel {num(total, 'indel'):,})")
    emit("\n== effect class (AGREE = snpEff gives a term that means the qmutfind label; CONVENTION = known difference of definition; DIFFER = other)")
    for c in ("SNV", "indel"):
        t = num(total, c)
        emit(f"  {c}: agree {num(agree, c):,} ({num(agree, c) / t:.3%}), convention {num(conv, c):,} ({num(conv, c) / t:.3%}), differ {num(differ, c):,} ({num(differ, c) / t:.3%}) of {t:,}")
    emit("\n  by qmutfind label: total, agree, convention, differ")
    for (c, label), t in sorted(total.items(), key=lambda kv: (kv[0][0], -kv[1])):
        emit(f"    {c:5s} {label:20s} {t:9,} {agree[(c, label)]:9,} {conv[(c, label)]:7,} {differ[(c, label)]:6,}")
    emit("\n== conventions seen")
    for r, n in why.most_common(): emit(f"    {n:7,}  {r}")
    emit("\n== DIFFER: qmutfind label -> snpEff term (the 25 most common)")
    for (a, b), n in first_term.most_common(25): emit(f"    {n:7,}  qmutfind {a:20s} snpEff {b}")
    emit("\n== protein change (variants that agree on class and where both give a change; HGVS.p turned into one-letter form)")
    emit("   exact comparison for missense, nonsense and in-frame changes; for start_lost, stop_lost and frameshift the residue and position only")
    for c in ("SNV", "indel"):
        t = num(aa_total, c); e = num(aa_equal, c)
        emit(f"  {c}: {e:,} of {t:,} equal = {e / max(t, 1):.4%}")
    for (c, label), t in sorted(aa_total.items()):
        emit(f"    {c:5s} {label:20s} {aa_equal[(c, label)]:9,} / {t:9,}")
    emit("\n== gene and transcript (variants that qmutfind puts in a gene)")
    for c in ("SNV", "indel"):
        t = gene_total[c]
        emit(f"  {c}: same transcript {gene_same[c]:,} of {t:,} ({gene_same[c] / max(t, 1):.4%}); snpEff names the same gene but another transcript or lists it later {gene_explained[c]:,}; other {t - gene_same[c] - gene_explained[c]:,}")
    open(sys.argv[3] + ".report.txt", "w").write("\n".join(L) + "\n")
    with open(sys.argv[3] + ".disagreements.tsv", "w") as fh:
        fh.write("chrom\tpos\tref\talt\tqmutfind\tsnpeff\tqmutfind_gene\tsnpeff_gene\tqmutfind_tx\tsnpeff_tx\tqmutfind_aa\tsnpeff_hgvs_p\tkind\n")
        for r in rows: fh.write("\t".join(map(str, r)) + "\n")
    with open(sys.argv[3] + ".aa_disagreements.tsv", "w") as fh:
        fh.write("chrom\tpos\tref\talt\tlabel\tqmutfind_aa\tsnpeff_aa\tsnpeff_hgvs_p\n")
        for r in aa_bad: fh.write("\t".join(map(str, r)) + "\n")
    with open(sys.argv[3] + ".gene_other.tsv", "w") as fh:
        fh.write("chrom\tpos\tref\talt\tlabel\tqmutfind_gene\tqmutfind_tx\tsnpeff_gene\tsnpeff_tx\tsnpeff_term\n")
        for r in gene_bad: fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()

"""Compare the effect labels of `qmutfind annotate` with `bcftools csq`.

Both tools annotate the same variants. The comparison is made on the canonical transcript that
qmutfind chose, because `bcftools csq` reports every transcript of a gene. Usage:

    python3 compare_annotation.py --mine annotated.vcf --csq csq.vcf [--examples 15]
"""
import argparse
import re
from collections import Counter

# My label -> the Sequence Ontology terms of bcftools csq that agree with it. A label that is not
# here (complex, large_indel_in_cds, unannotated) has no counterpart and is not compared.
LABEL_TO_TERMS = {
    "synonymous": {"synonymous", "stop_retained", "start_retained"},
    "missense": {"missense"},
    "nonsense": {"stop_gained"},
    "start_lost": {"start_lost"},
    "stop_lost": {"stop_lost"},
    "frameshift": {"frameshift"},
    "inframe_insertion": {"inframe_insertion", "inframe_altering"},
    "inframe_deletion": {"inframe_deletion", "inframe_altering"},
    "splice_site": {"splice_donor", "splice_acceptor"},
    "intron_derived": {"intron"},
    "5UTR_derived": {"5_prime_utr"},
    "3UTR_derived": {"3_prime_utr"},
    "non_coding": {"non_coding"},
}
PROTEIN_SNP_LABELS = {"synonymous", "missense", "nonsense", "start_lost", "stop_lost"}


def parse_bcsq(info):
    """The entries of the BCSQ field of one INFO string. An entry that starts with `@` only points
    at another record and is dropped."""
    for part in info.split(";"):
        if part.startswith("BCSQ="):
            entries = []
            for raw in part[len("BCSQ="):].split(","):
                if raw.startswith("@"):
                    continue
                fields = raw.split("|") + [""] * 7
                entries.append({
                    "terms": [t for t in fields[0].lstrip("*").split("&") if t],
                    "gene": fields[1],
                    "transcript": fields[2],
                    "biotype": fields[3],
                    "strand": fields[4],
                    "aa": fields[5],
                    "dna": fields[6],
                })
            return entries
    return []


def parse_csq_aa(text):
    """(reference residue, position, alternate residue) from a bcftools csq amino-acid field, or None."""
    m = re.fullmatch(r"(\d+)([A-Z*])>(\d+)([A-Z*])", text or "")
    if m and m.group(1) == m.group(3):
        return m.group(2), int(m.group(1)), m.group(4)
    m = re.fullmatch(r"(\d+)([A-Z*])", text or "")  # synonymous: one residue
    if m:
        return m.group(2), int(m.group(1)), m.group(2)
    return None


def parse_my_aa(text):
    """(reference residue, position, alternate residue) from a qmutfind AA tag, or None."""
    m = re.fullmatch(r"([A-Z*])(\d+)([A-Z*])", text or "")
    return (m.group(1), int(m.group(2)), m.group(3)) if m else None


def effect_matches(label, terms):
    return bool(LABEL_TO_TERMS.get(label, set()) & set(terms))


def compare_annotations(mine, csq, max_examples=50, gene_aliases=None):
    """`mine` maps a variant key to {effect, transcript, aa}. `csq` maps a key to the parsed BCSQ entries."""
    gene_aliases = gene_aliases or {}
    r = Counter()
    confusion = Counter()
    examples = []

    def example(kind, key, mine_rec, entry):
        if len(examples) < max_examples:
            examples.append({
                "kind": kind, "key": key, "mine": mine_rec["effect"], "mine_aa": mine_rec["aa"],
                "csq": "&".join(entry["terms"]) if entry else None, "csq_aa": entry["aa"] if entry else None,
            })

    for key, rec in mine.items():
        r["total"] += 1
        entries = csq.get(key, [])
        label, tx = rec["effect"], rec["transcript"]
        if label == "intergenic":
            if not entries:
                r["both_silent"] += 1
                confusion[(label, "none")] += 1
            else:
                r["mine_intergenic_csq_genic"] += 1
                confusion[(label, entries[0]["terms"][0] if entries[0]["terms"] else "none")] += 1
                example("mine_intergenic_csq_genic", key, rec, entries[0])
            continue
        if not entries:
            r["csq_silent_mine_genic"] += 1
            confusion[(label, "none")] += 1
            example("csq_silent", key, rec, None)
            continue
        # bcftools csq leaves the transcript empty for introns, so those entries are matched by gene.
        names = {rec.get("gene")} | gene_aliases.get(rec.get("gene"), set())
        on_tx = [e for e in entries if e["transcript"] == tx or (not e["transcript"] and e.get("gene") and e["gene"] in names)]
        if not on_tx:
            r["csq_other_transcript"] += 1
            confusion[(label, "other_transcript")] += 1
            continue
        terms = {t for e in on_tx for t in e["terms"]}
        confusion[(label, on_tx[0]["terms"][0] if on_tx[0]["terms"] else "none")] += 1
        if label not in LABEL_TO_TERMS:
            r["not_comparable"] += 1
        elif label in ("5UTR_derived", "3UTR_derived") and terms == {"intron"}:
            # The GFF has no UTR features, so csq cannot tell a UTR from an intron.
            r["utr_vs_csq_intron"] += 1
        elif effect_matches(label, terms):
            r["agree"] += 1
        else:
            r["disagree"] += 1
            example("disagree", key, rec, on_tx[0])
        if label in PROTEIN_SNP_LABELS:
            mine_aa = parse_my_aa(rec["aa"])
            csq_aa = next((a for a in (parse_csq_aa(e["aa"]) for e in on_tx) if a), None)
            if mine_aa and csq_aa:
                r["aa_compared"] += 1
                if mine_aa == csq_aa:
                    r["aa_agree"] += 1
                else:
                    r["aa_disagree"] += 1
                    example("aa_disagree", key, rec, on_tx[0])
    out = dict(r)
    for name in ("total", "agree", "disagree", "csq_other_transcript", "both_silent", "csq_silent_mine_genic",
                 "mine_intergenic_csq_genic", "not_comparable", "utr_vs_csq_intron", "aa_compared", "aa_agree", "aa_disagree"):
        out.setdefault(name, 0)
    out["confusion"] = dict(confusion)
    out["examples"] = examples
    return out


def _info(info_field, key):
    for part in info_field.split(";"):
        if part.startswith(key + "="):
            return part[len(key) + 1:]
    return None


def load_mine(path):
    mine = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            effect = (_info(f[7], "EFFECT") or "unannotated").split(",")[0]
            aa = (_info(f[7], "AA") or "").split(",")[0]
            mine[(f[0], int(f[1]), f[3], f[4])] = {
                "effect": effect, "transcript": _info(f[7], "TRANSCRIPT"), "gene": _info(f[7], "GENE"), "aa": None if aa in ("", ".") else aa,
            }
    return mine


def load_csq(path):
    csq = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            csq[(f[0], int(f[1]), f[3], f[4])] = parse_bcsq(f[7])
    return csq


def load_gene_aliases(gff_path):
    """locus_tag -> {Name, gene} for the gene lines of a GFF3 (plain or gzip)."""
    import gzip
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    aliases = {}
    with opener(gff_path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                attrs = dict(kv.split("=", 1) for kv in f[8].strip().split(";") if "=" in kv)
                tag = attrs.get("locus_tag")
                if tag:
                    aliases[tag] = {v for k, v in attrs.items() if k in ("Name", "gene") and v != tag}
    return aliases


def main(argv=None):
    ap = argparse.ArgumentParser(description="Compare qmutfind annotate with bcftools csq.")
    ap.add_argument("--mine", required=True)
    ap.add_argument("--csq", required=True)
    ap.add_argument("--examples", type=int, default=15)
    ap.add_argument("--gff", help="a GFF3 for the gene Name aliases that bcftools csq may write")
    args = ap.parse_args(argv)
    aliases = load_gene_aliases(args.gff) if args.gff else None
    r = compare_annotations(load_mine(args.mine), load_csq(args.csq), max_examples=max(args.examples * 4, 50), gene_aliases=aliases)
    for name in ("total", "agree", "disagree", "csq_other_transcript", "both_silent", "csq_silent_mine_genic",
                 "mine_intergenic_csq_genic", "not_comparable", "utr_vs_csq_intron", "aa_compared", "aa_agree", "aa_disagree"):
        print(f"{name}\t{r[name]}")
    print("\nmine\tcsq\tcount")
    for (a, b), n in sorted(r["confusion"].items(), key=lambda x: -x[1]):
        print(f"{a}\t{b}\t{n}")
    shown = Counter()
    print("\nkind\tvariant\tmine\tmine_aa\tcsq\tcsq_aa")
    for e in r["examples"]:
        if shown[e["kind"]] < args.examples:
            shown[e["kind"]] += 1
            k = e["key"]
            print(f"{e['kind']}\t{k[0]}:{k[1]} {k[2][:12]}>{k[3][:12]}\t{e['mine']}\t{e['mine_aa']}\t{e['csq']}\t{e['csq_aa']}")


if __name__ == "__main__":
    main()

"""Design known mutations for a benchmark truth set.

Positions are 1-based and in reference coordinates, as in a VCF.
An indel record carries one anchor base, as the VCF format requires.
"""
import argparse
import gzip
import random
import subprocess
from contextlib import contextmanager
from dataclasses import dataclass

# ASSUMED defaults. They are not derived from data. See PLAN.md section 11.
SMALL_INDEL_LEN = (1, 10)  # the small class is 1 to 49 bp
LARGE_INDEL_LEN = (50, 500)
END_MARGIN = 1000  # keep variants away from contig ends
MIN_SPACING = 200  # bp between the end of one variant and the start of the next


@dataclass(frozen=True)
class Variant:
    contig: str
    pos: int
    ref: str
    alt: str
    klass: str  # SNP, SMALL_INDEL or LARGE_INDEL


def _pick_contig(reference, rng):
    names = list(reference)
    return rng.choices(names, weights=[len(reference[n]) for n in names])[0]


def _make_variant(reference, klass, rng, site=None):
    """One random variant of class `klass`. `site` is (contig, pos) to place it at a chosen position."""
    if site is None:
        contig = _pick_contig(reference, rng)
        pos = rng.randint(END_MARGIN, len(reference[contig]) - END_MARGIN - LARGE_INDEL_LEN[1])
    else:
        contig, pos = site
    seq = reference[contig]
    if klass == "SNP":
        ref = seq[pos - 1]
        alt = rng.choice([b for b in "ACGT" if b != ref])
        return Variant(contig, pos, ref, alt, klass)
    low, high = SMALL_INDEL_LEN if klass == "SMALL_INDEL" else LARGE_INDEL_LEN
    length = rng.randint(low, high)
    anchor = seq[pos - 1]
    if rng.random() < 0.5:  # deletion
        return Variant(contig, pos, seq[pos - 1 : pos + length], anchor, klass)
    inserted = "".join(rng.choice("ACGT") for _ in range(length))
    return Variant(contig, pos, anchor, anchor + inserted, klass)


def apply_variants(reference, variants):
    """Return the mutated sequences. Apply from the highest position down so that
    the coordinates of the variants not yet applied stay valid."""
    mutated = {}
    for contig, seq in reference.items():
        here = sorted((v for v in variants if v.contig == contig), key=lambda v: v.pos, reverse=True)
        for v in here:
            start = v.pos - 1
            seq = seq[:start] + v.alt + seq[start + len(v.ref):]
        mutated[contig] = seq
    return mutated


@contextmanager
def _open_text(path):
    """Open a plain, .gz or .zst text file for reading. .zst uses `zstd -dc`."""
    path = str(path)
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            yield fh
    elif path.endswith(".zst"):
        proc = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True)
        try:
            yield proc.stdout
        finally:
            proc.stdout.close()
            if proc.wait() != 0:
                raise RuntimeError(f"zstd failed on {path}")
    else:
        with open(path) as fh:
            yield fh


def read_fasta(path):
    """Return {contig name: uppercase sequence}. The name is the first word of the header."""
    seqs, name, parts = {}, None, []
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts).upper()
                name, parts = line[1:].split()[0], []
            elif name is not None:
                parts.append(line.strip())
    if name is not None:
        seqs[name] = "".join(parts).upper()
    return seqs


def write_fasta(seqs, out, width=80):
    """Write {name: sequence} as FASTA to the text handle `out`."""
    for name, seq in seqs.items():
        out.write(f">{name}\n")
        for i in range(0, len(seq), width):
            out.write(seq[i : i + width] + "\n")


def write_truth_vcf(variants, reference, out):
    """Write the truth set as VCF 4.2 to the text handle `out`.
    Records are sorted by contig order in `reference`, then by position."""
    out.write("##fileformat=VCFv4.2\n")
    for contig, seq in reference.items():
        out.write(f"##contig=<ID={contig},length={len(seq)}>\n")
    out.write('##INFO=<ID=CLASS,Number=1,Type=String,Description="SNP, SMALL_INDEL or LARGE_INDEL">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    order = {name: i for i, name in enumerate(reference)}
    for v in sorted(variants, key=lambda v: (order[v.contig], v.pos)):
        out.write(f"{v.contig}\t{v.pos}\t.\t{v.ref}\t{v.alt}\t.\tPASS\tCLASS={v.klass}\n")


def _usable(candidate, reference):
    """True when the reference footprint holds only A, C, G and T."""
    start = candidate.pos - 1
    footprint = reference[candidate.contig][start : start + len(candidate.ref)]
    return set(footprint) <= set("ACGT")


def _too_close(candidate, taken, spacing=MIN_SPACING):
    start, end = candidate.pos, candidate.pos + len(candidate.ref) - 1
    for s, e in taken.get(candidate.contig, ()):
        if not (start - e >= spacing or s - end >= spacing):
            return True
    return False


def _homopolymer_indel(reference, rng):
    """A one-base indel in a run of four or more identical bases: one more or one fewer copy of the base.
    The VCF anchor is the base before the run. Returns None when no run is found near the random start."""
    contig = _pick_contig(reference, rng)
    seq = reference[contig]
    start = rng.randint(END_MARGIN, len(seq) - END_MARGIN - 5_000)
    for i in range(start, start + 5_000):
        base = seq[i]
        if base in "ACGT" and seq[i : i + 4] == base * 4 and seq[i - 1] != base:
            anchor = seq[i - 1]
            if rng.random() < 0.5:  # delete one copy of the base
                return Variant(contig, i, anchor + base, anchor, "SMALL_INDEL")
            return Variant(contig, i, anchor, anchor + base, "SMALL_INDEL")  # add one copy
    return None


def _cluster_site(taken, rng):
    """A (contig, pos) 6 to 40 bp from a variant that is already placed."""
    contig = rng.choice(sorted(taken))
    start, end = rng.choice(taken[contig])
    offset = rng.randint(6, 40)
    return contig, (end + offset if rng.random() < 0.5 else max(1, start - offset))


def design_variants(reference, n_snp, n_small, n_large, seed, max_attempts=1000,
                    min_spacing=MIN_SPACING, homopolymer_fraction=0.0, cluster_fraction=0.0):
    """Design the variants. `homopolymer_fraction` of the small indels are one-base indels in
    homopolymer runs. `cluster_fraction` of the other variants are placed 6 to 40 bp from an earlier
    one. Both make the truth set harder than uniform random positions."""
    rng = random.Random(seed)
    n_homopolymer = round(homopolymer_fraction * n_small)
    plan = (
        [("SNP", False)] * n_snp
        + [("SMALL_INDEL", True)] * n_homopolymer
        + [("SMALL_INDEL", False)] * (n_small - n_homopolymer)
        + [("LARGE_INDEL", False)] * n_large
    )
    taken = {}
    variants = []
    for klass, homopolymer in plan:
        placed = None
        for _ in range(max_attempts):
            clustered = (not homopolymer) and bool(variants) and rng.random() < cluster_fraction
            if homopolymer:
                cand = _homopolymer_indel(reference, rng)
                if cand is None:
                    continue
            else:
                site = _cluster_site(taken, rng) if clustered else None
                cand = _make_variant(reference, klass, rng, site)
            if cand.pos < END_MARGIN or cand.pos + len(cand.ref) > len(reference[cand.contig]) - END_MARGIN:
                continue
            if _usable(cand, reference) and not _too_close(cand, taken, 5 if clustered else min_spacing):
                placed = cand
                break
        if placed is None:
            raise ValueError(f"could not place a {klass} after {max_attempts} attempts (spacing {min_spacing})")
        taken.setdefault(placed.contig, []).append((placed.pos, placed.pos + len(placed.ref) - 1))
        variants.append(placed)
    return variants


def main(argv=None):
    ap = argparse.ArgumentParser(description="Design a truth set and write the mutated genome.")
    ap.add_argument("--reference", required=True, help="reference FASTA (plain, .gz or .zst)")
    ap.add_argument("--out-prefix", required=True, help="writes <prefix>.mutated.fa and <prefix>.truth.vcf")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--snps", type=int, default=0)
    ap.add_argument("--small-indels", type=int, default=0)
    ap.add_argument("--large-indels", type=int, default=0)
    ap.add_argument("--min-spacing", type=int, default=MIN_SPACING, help="bp between variants (default 200)")
    ap.add_argument("--homopolymer-fraction", type=float, default=0.0, help="share of small indels placed in homopolymer runs")
    ap.add_argument("--cluster-fraction", type=float, default=0.0, help="share of variants placed close to an earlier one")
    args = ap.parse_args(argv)

    reference = read_fasta(args.reference)
    variants = design_variants(
        reference, args.snps, args.small_indels, args.large_indels, args.seed,
        min_spacing=args.min_spacing, homopolymer_fraction=args.homopolymer_fraction,
        cluster_fraction=args.cluster_fraction,
    )
    with open(f"{args.out_prefix}.mutated.fa", "w") as fh:
        write_fasta(apply_variants(reference, variants), fh)
    with open(f"{args.out_prefix}.truth.vcf", "w") as fh:
        write_truth_vcf(variants, reference, fh)


if __name__ == "__main__":
    main()

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


def _make_variant(reference, klass, rng):
    contig = _pick_contig(reference, rng)
    seq = reference[contig]
    pos = rng.randint(END_MARGIN, len(seq) - END_MARGIN - LARGE_INDEL_LEN[1])
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


def _too_close(candidate, taken):
    start, end = candidate.pos, candidate.pos + len(candidate.ref) - 1
    for s, e in taken.get(candidate.contig, ()):
        if not (start - e >= MIN_SPACING or s - end >= MIN_SPACING):
            return True
    return False


def design_variants(reference, n_snp, n_small, n_large, seed, max_attempts=1000):
    rng = random.Random(seed)
    klasses = ["SNP"] * n_snp + ["SMALL_INDEL"] * n_small + ["LARGE_INDEL"] * n_large
    taken = {}
    variants = []
    for klass in klasses:
        for _ in range(max_attempts):
            cand = _make_variant(reference, klass, rng)
            if _usable(cand, reference) and not _too_close(cand, taken):
                break
        else:
            raise ValueError(
                f"could not place a {klass} with spacing {MIN_SPACING} after {max_attempts} attempts"
            )
        taken.setdefault(cand.contig, []).append((cand.pos, cand.pos + len(cand.ref) - 1))
        variants.append(cand)
    return variants


def main(argv=None):
    ap = argparse.ArgumentParser(description="Design a truth set and write the mutated genome.")
    ap.add_argument("--reference", required=True, help="reference FASTA (plain, .gz or .zst)")
    ap.add_argument("--out-prefix", required=True, help="writes <prefix>.mutated.fa and <prefix>.truth.vcf")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--snps", type=int, default=0)
    ap.add_argument("--small-indels", type=int, default=0)
    ap.add_argument("--large-indels", type=int, default=0)
    args = ap.parse_args(argv)

    reference = read_fasta(args.reference)
    variants = design_variants(
        reference, args.snps, args.small_indels, args.large_indels, args.seed
    )
    with open(f"{args.out_prefix}.mutated.fa", "w") as fh:
        write_fasta(apply_variants(reference, variants), fh)
    with open(f"{args.out_prefix}.truth.vcf", "w") as fh:
        write_truth_vcf(variants, reference, fh)


if __name__ == "__main__":
    main()

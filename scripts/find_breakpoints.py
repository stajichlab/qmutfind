#!/usr/bin/env python3
"""Find assembly-vs-reference breakpoints near a locus with minimap2.

For each assembly FASTA (plain, .gz or .zst):
  1. Run `minimap2 -x asm5 --cs` against the reference (skipped if the PAF exists).
  2. Keep alignments on the target scaffold that overlap LOCUS +/- FLANK.
  3. Order same-strand blocks by reference position and report each gap
     between consecutive blocks: reference breakpoints, query gap size,
     reference gap size, and a call (insertion / deletion / replacement).
  4. Write the query sequence inside each gap to a FASTA file.

Usage:
  find_breakpoints.py --ref A1163.fasta.gz --gff A1163.gff.gz --gene AFUB_101460 \
      --outdir results/breakpoints results/sort/*.sorted.fasta.gz
  find_breakpoints.py --ref ref.fa --locus scf_000011_A_fumigatus_A1163:282234-285136 ...

minimap2 must be on PATH (module load minimap2).
"""
import argparse
import gzip
import os
import re
import subprocess
import sys

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    if path.endswith(".zst"):
        p = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True)
        return p.stdout
    return open(path)


def sample_name(path):
    b = os.path.basename(path)
    return re.sub(r"\.(sorted\.)?(fasta|fa|fna)(\.gz|\.zst)?$", "", b)


def locus_from_gff(gff, gene):
    with open_text(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "protein_coding_gene"):
                continue
            if re.search(r"(^|;)ID=%s(;|$)" % re.escape(gene), f[8]):
                return f[0], int(f[3]), int(f[4])
    sys.exit("gene %s not found in %s" % (gene, gff))


def read_fasta(path, wanted):
    seqs, name, buf = {}, None, []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name in wanted:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name in wanted:
                buf.append(line.strip())
    if name in wanted:
        seqs[name] = "".join(buf)
    return seqs


def run_minimap2(ref, asm, paf, threads, preset):
    if os.path.exists(paf) and os.path.getsize(paf) > 0:
        return
    cmd = ["minimap2", "-x", preset, "-t", str(threads), "--cs", ref, asm]
    with open(paf + ".tmp", "w") as out:
        r = subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL)
    if r.returncode != 0:
        sys.exit("minimap2 failed for %s" % asm)
    os.replace(paf + ".tmp", paf)


def blocks_in_window(paf, scaf, lo, hi, min_aln):
    by_query = {}
    with open(paf) as fh:
        for line in fh:
            f = line.split("\t")
            if f[5] != scaf or int(f[10]) < min_aln:
                continue
            ts, te = int(f[7]), int(f[8])
            if te < lo or ts > hi:
                continue
            by_query.setdefault(f[0], []).append(
                dict(q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]),
                     strand=f[4], ts=ts, te=te, alen=int(f[10]), mapq=int(f[11])))
    return by_query


def call_gaps(by_query, gs, ge):
    rows = []
    for q, bl in by_query.items():
        for strand in ("+", "-"):
            sb = sorted([b for b in bl if b["strand"] == strand], key=lambda b: b["ts"])
            for a, b in zip(sb, sb[1:]):
                rgap = b["ts"] - a["te"]
                qgap = (b["qs"] - a["qe"]) if strand == "+" else (a["qs"] - b["qe"])
                if qgap < 0 or rgap < 0:  # overlapping blocks: not a clean gap
                    continue
                if qgap == 0 and rgap == 0:
                    continue
                if rgap == 0:
                    kind = "insertion"
                elif qgap == 0:
                    kind = "deletion"
                else:
                    kind = "replacement"
                q_lo = (a["qe"] if strand == "+" else b["qe"])
                q_hi = (b["qs"] if strand == "+" else a["qs"])
                rows.append(dict(
                    contig=q, strand=strand, ref_left_end=a["te"], ref_right_start=b["ts"] + 1,
                    ref_gap=rgap, query_gap=qgap, q_lo=q_lo, q_hi=q_hi, kind=kind,
                    in_locus=(a["te"] < ge and b["ts"] + 1 > gs),
                    mapq_min=min(a["mapq"], b["mapq"])))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("assemblies", nargs="+")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--locus", help="scaffold:start-end (1-based)")
    ap.add_argument("--gff")
    ap.add_argument("--gene")
    ap.add_argument("--flank", type=int, default=50000)
    ap.add_argument("--min-aln", type=int, default=1000)
    ap.add_argument("--min-gap", type=int, default=20, help="ignore gaps smaller than this on both axes")
    ap.add_argument("--preset", default="asm5")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--outdir", default="breakpoints")
    a = ap.parse_args()

    if a.locus:
        m = re.match(r"(.+):(\d+)-(\d+)$", a.locus)
        if not m:
            sys.exit("--locus must look like scaffold:start-end")
        scaf, gs, ge = m.group(1), int(m.group(2)), int(m.group(3))
    elif a.gff and a.gene:
        scaf, gs, ge = locus_from_gff(a.gff, a.gene)
    else:
        sys.exit("give --locus, or --gff and --gene")

    os.makedirs(a.outdir, exist_ok=True)
    summary = os.path.join(a.outdir, "breakpoint_summary.tsv")
    cols = ["sample", "contig", "strand", "ref_scaffold", "ref_left_end", "ref_right_start",
            "ref_gap_bp", "query_gap_bp", "call", "overlaps_locus", "min_mapq", "insert_fasta_id"]
    with open(summary, "w") as sf, open(os.path.join(a.outdir, "inserted_sequences.fasta"), "w") as ff:
        sf.write("\t".join(cols) + "\n")
        for asm in a.assemblies:
            s = sample_name(asm)
            paf = os.path.join(a.outdir, s + ".vs_ref.paf")
            run_minimap2(a.ref, asm, paf, a.threads, a.preset)
            by_q = blocks_in_window(paf, scaf, gs - a.flank, ge + a.flank, a.min_aln)
            rows = [r for r in call_gaps(by_q, gs, ge)
                    if max(r["ref_gap"], r["query_gap"]) >= a.min_gap]
            if not by_q:
                print("%s: no alignments within %s:%d-%d" % (s, scaf, gs - a.flank, ge + a.flank))
                continue
            seqs = read_fasta(asm, {r["contig"] for r in rows})
            for r in rows:
                fid = ""
                if r["query_gap"] > 0 and r["contig"] in seqs:
                    seq = seqs[r["contig"]][r["q_lo"]:r["q_hi"]]
                    if r["strand"] == "-":
                        seq = seq.translate(COMP)[::-1]
                    fid = "%s|%s:%d-%d|ref_%d_%d" % (s, r["contig"], r["q_lo"] + 1, r["q_hi"],
                                                    r["ref_left_end"], r["ref_right_start"])
                    ff.write(">%s\n%s\n" % (fid, seq))
                sf.write("\t".join(map(str, [
                    s, r["contig"], r["strand"], scaf, r["ref_left_end"], r["ref_right_start"],
                    r["ref_gap"], r["query_gap"], r["kind"], "yes" if r["in_locus"] else "no",
                    r["mapq_min"], fid])) + "\n")
                if r["in_locus"]:
                    print("%s: %s of %d bp (ref lost %d bp) between %s:%d and %d on %s" % (
                        s, r["kind"], r["query_gap"], r["ref_gap"], scaf, r["ref_left_end"],
                        r["ref_right_start"], r["contig"]))
            if not any(r["in_locus"] for r in rows):
                print("%s: alignments present, no gap overlapping %s:%d-%d" % (s, scaf, gs, ge))
    print("wrote %s" % summary)


if __name__ == "__main__":
    main()

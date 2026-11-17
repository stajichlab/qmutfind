"""scripts/find_breakpoints.py: it finds where an assembly differs from the reference near a locus, with minimap2, and writes the table
that `qmutfind locus` and confirm_allele.html read. The assembly is simulated, so the right answer is known."""
import csv
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "find_breakpoints.py"
pytestmark = pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 is not installed (run through pixi)")

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def random_dna(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def fasta(path, name, seq):
    path.write_text(f">{name}\n{seq}\n")
    return path


def run(tmp_path, assembly, extra, name="sample"):
    out = tmp_path / "out"
    cmd = [sys.executable, "-I", str(SCRIPT), "--ref", str(tmp_path / "ref.fa"), "--outdir", str(out), *extra, str(assembly)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader((out / "breakpoint_summary.tsv").open(), delimiter="\t"))
    inserted = {}
    key = None
    for line in (out / "inserted_sequences.fasta").read_text().splitlines():
        if line.startswith(">"):
            key = line[1:]
            inserted[key] = ""
        else:
            inserted[key] += line
    return result.stdout, rows, inserted


@pytest.fixture
def world(tmp_path):
    genome = random_dna(60_000, 1)
    insert = random_dna(1_500, 2)
    fasta(tmp_path / "ref.fa", "chr", genome)
    # 20,001-22,500 (2,500 bases) are replaced by 1,500 other bases.
    assembly = genome[:20_000] + insert + genome[22_500:]
    (tmp_path / "genes.gff3").write_text("##gff-version 3\nchr\ts\tgene\t19001\t23000\t.\t+\t.\tID=G1\n")
    return tmp_path, genome, insert, assembly


def near(value, expected, slack=8):
    # A few bases at a junction can equal the reference by chance, which moves the breakpoint by that many bases.
    return abs(int(value) - expected) <= slack


def test_a_replacement_near_a_locus_is_found_with_its_flanks_sizes_and_inserted_sequence(world):
    tmp_path, genome, insert, assembly = world
    fasta(tmp_path / "mut.fa", "ctg1", assembly)

    stdout, rows, inserted = run(tmp_path, tmp_path / "mut.fa", ["--locus", "chr:19001-23000"])

    hits = [r for r in rows if r["overlaps_locus"] == "yes"]
    assert len(hits) == 1, rows
    row = hits[0]
    assert row["call"] == "replacement" and row["ref_scaffold"] == "chr" and row["contig"] == "ctg1"
    assert near(row["ref_left_end"], 20_000) and near(row["ref_right_start"], 22_501)
    assert near(row["query_gap_bp"], 1_500) and near(row["ref_gap_bp"], 2_500)
    assert int(row["ref_gap_bp"]) == int(row["ref_right_start"]) - int(row["ref_left_end"]) - 1
    seq = inserted[row["insert_fasta_id"]]
    assert len(seq) == int(row["query_gap_bp"])
    assert seq[3:-3] in insert, "the inserted sequence is the simulated one (a few junction bases can equal the reference)"
    assert "replacement of" in stdout and "between chr:" in stdout


def test_an_assembly_on_the_other_strand_gives_the_same_flanks_and_the_inserted_sequence_in_reference_orientation(world):
    tmp_path, genome, insert, assembly = world
    fasta(tmp_path / "rev.fa", "ctg1", revcomp(assembly))

    _, rows, inserted = run(tmp_path, tmp_path / "rev.fa", ["--locus", "chr:19001-23000"])

    row = [r for r in rows if r["overlaps_locus"] == "yes"][0]
    assert row["strand"] == "-"
    assert near(row["ref_left_end"], 20_000) and near(row["ref_right_start"], 22_501)
    seq = inserted[row["insert_fasta_id"]]
    assert seq[5:-5] in insert, "the sequence is written as it reads on the reference strand"


def test_a_gene_id_from_a_gff_gives_the_locus_and_an_unchanged_assembly_has_no_gap(world):
    tmp_path, genome, insert, assembly = world
    fasta(tmp_path / "wt.fa", "ctg1", genome)
    fasta(tmp_path / "mut.fa", "ctg1", assembly)

    stdout, rows, _ = run(tmp_path, tmp_path / "wt.fa", ["--gff", str(tmp_path / "genes.gff3"), "--gene", "G1"])
    assert rows == [] and "no gap overlapping chr:19001-23000" in stdout

    _, rows, _ = run(tmp_path, tmp_path / "mut.fa", ["--gff", str(tmp_path / "genes.gff3"), "--gene", "G1"])
    assert [r["overlaps_locus"] for r in rows] == ["yes"]


def test_a_missing_locus_or_gene_is_an_error_that_says_what_to_give(world):
    tmp_path, *_ = world
    fasta(tmp_path / "mut.fa", "ctg1", "ACGT" * 100)
    result = subprocess.run([sys.executable, "-I", str(SCRIPT), "--ref", str(tmp_path / "ref.fa"), str(tmp_path / "mut.fa")], capture_output=True, text=True)
    assert result.returncode != 0 and "--locus" in result.stderr

import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "snpeff_validation"))
from compare_snpeff import ACCEPT, convention, one_letter, ref_and_position  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "snpeff_validation" / "compare_snpeff.py"


def test_hgvs_protein_changes_become_one_letter_form():
    assert one_letter("p.Val29Ile") == "V29I"
    assert one_letter("p.Trp8Ter") == "W8*"
    assert one_letter("p.Glu3Aspfs*4") == "E3Dfs*4"
    assert one_letter("p.Met1?") == "M1?"
    assert one_letter("p.Ter244Argext*?") == "*244Rext*?"
    assert one_letter("p.Leu29=") == "L29="
    assert one_letter("c.85G>A") == ""
    assert one_letter("") == ""


def test_the_residue_and_position_are_read_from_either_notation():
    assert ref_and_position("M1V") == ("M", "1")
    assert ref_and_position("*244R") == ("*", "244")
    assert ref_and_position("E3Dfs*4") == ("E", "3")
    assert ref_and_position("?471Afs*?") is None
    assert ref_and_position("") is None


def test_every_label_that_annotate_writes_has_a_set_of_snpeff_terms():
    labels = {"synonymous", "missense", "nonsense", "start_lost", "stop_lost", "frameshift", "inframe_insertion", "inframe_deletion", "splice_site",
              "5UTR_derived", "3UTR_derived", "intron_derived", "non_coding", "complex", "large_indel_in_cds", "intergenic"}
    assert labels <= set(ACCEPT)
    assert ACCEPT["missense"] == {"missense_variant"} and "stop_gained" in ACCEPT["nonsense"]


@pytest.mark.parametrize("label, terms, expected", [
    ("start_lost", {"initiator_codon_variant"}, "initiator_codon_variant"),
    ("5UTR_derived", {"splice_region_variant"}, "splice_region_variant"),
    ("intron_derived", {"splice_region_variant"}, "splice_region_variant"),
    ("non_coding", {"intragenic_variant"}, "intragenic_variant"),
    ("missense", {"frameshift_variant"}, None),
    ("splice_site", {"splice_region_variant", "intron_variant"}, None),  # not only splice_region_variant: a real difference
    ("intron_derived", {"intragenic_variant"}, None),
])
def test_only_the_three_known_conventions_are_called_conventions(label, terms, expected):
    result = convention(label, terms)
    assert (result is None) if expected is None else (expected in result)


def _vcf(path, rows, ann=False):
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for chrom, pos, ref, alt, info in rows:
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}")
    path.write_text("\n".join(lines) + "\n")


def test_the_comparison_counts_agree_convention_and_differ_on_a_small_pair_of_files(tmp_path):
    q, s = tmp_path / "q.vcf", tmp_path / "s.vcf"
    _vcf(q, [
        ("c", 10, "A", "C", "EFFECT=missense;GENE=G1;TRANSCRIPT=tx1;AA=V29I"),
        ("c", 20, "A", "C", "EFFECT=missense;GENE=G1;TRANSCRIPT=tx1;AA=V30I"),
        ("c", 30, "A", "C", "EFFECT=start_lost;GENE=G1;TRANSCRIPT=tx1;AA=M1L"),
        ("c", 40, "A", "C", "EFFECT=intergenic"),
        ("c", 50, "a", "ac", "EFFECT=intergenic"),  # lower case in the reference allele: compared as upper case
    ])
    _vcf(s, [
        ("c", 10, "A", "C", "ANN=C|missense_variant|MODERATE|G1|G1|transcript|tx1|protein_coding|1/1|c.85A>C|p.Val29Ile|||||"),
        ("c", 20, "A", "C", "ANN=C|stop_gained|HIGH|G1|G1|transcript|tx1|protein_coding|1/1|c.86A>C|p.Val30Ter|||||"),
        ("c", 30, "A", "C", "ANN=C|initiator_codon_variant|LOW|G1|G1|transcript|tx1|protein_coding|1/1|c.1A>C|p.Met1?|||||"),
        ("c", 40, "A", "C", "ANN=C|intergenic_region|MODIFIER|G1-G2|G1-G2|intergenic_region|G1-G2||||||||"),
        ("c", 50, "A", "AC", "ANN=AC|intergenic_region|MODIFIER|G1-G2|G1-G2|intergenic_region|G1-G2||||||||"),
    ])
    out = subprocess.run([sys.executable, str(SCRIPT), str(q), str(s), str(tmp_path / "out")], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    report = (tmp_path / "out.report.txt").read_text()
    assert "variants compared: 5" in report
    assert "SNV: agree 2 (50.000%), convention 1 (25.000%), differ 1 (25.000%) of 4" in report
    assert "indel: agree 1 (100.000%)" in report
    assert "start_lost / snpEff initiator_codon_variant" in report
    rows = (tmp_path / "out.disagreements.tsv").read_text().splitlines()
    assert len(rows) == 3 and rows[1].endswith("differ") and rows[2].endswith("convention"), rows  # in the order of the positions: 20, then 30

import gzip
import math

from compare_calls import compare_variants, load_svlens, load_vcf, main

# A variant is keyed by (contig, pos, ref, alt). The value is its class.
TRUTH = {
    ("c", 10, "A", "G"): "SNP",
    ("c", 20, "A", "C"): "SNP",
    ("c", 30, "A", "T"): "SNP",
    ("c", 40, "AT", "A"): "SMALL_INDEL",
    ("c", 5000, "C" + "A" * 50, "C"): "LARGE_INDEL",
}
CALLS = {
    ("c", 10, "A", "G"): "SNP",          # true positive
    ("c", 20, "A", "C"): "SNP",          # true positive
    ("c", 99, "A", "G"): "SNP",          # false positive
    ("c", 40, "AT", "A"): "SMALL_INDEL",  # true positive
    ("c", 41, "T", "TG"): "SMALL_INDEL",  # false positive
}


def test_counts_per_class_are_true_positives_false_negatives_and_false_positives():
    stats = compare_variants(TRUTH, CALLS)

    assert stats["SNP"]["truth"] == 3
    assert stats["SNP"]["tp"] == 2
    assert stats["SNP"]["fn"] == 1
    assert stats["SNP"]["calls"] == 3
    assert stats["SNP"]["fp"] == 1
    assert stats["SMALL_INDEL"]["tp"] == 1 and stats["SMALL_INDEL"]["fp"] == 1
    assert stats["LARGE_INDEL"]["truth"] == 1 and stats["LARGE_INDEL"]["tp"] == 0


def test_sensitivity_and_precision_are_ratios_and_nan_when_the_denominator_is_zero():
    stats = compare_variants(TRUTH, CALLS)

    assert math.isclose(stats["SNP"]["sensitivity"], 2 / 3)
    assert math.isclose(stats["SNP"]["precision"], 2 / 3)
    assert math.isclose(stats["SMALL_INDEL"]["sensitivity"], 1.0)
    assert math.isclose(stats["SMALL_INDEL"]["precision"], 0.5)
    # The large indel class has truth but no calls: sensitivity 0, precision undefined.
    assert stats["LARGE_INDEL"]["sensitivity"] == 0.0
    assert math.isnan(stats["LARGE_INDEL"]["precision"])


def test_a_call_with_the_right_position_but_a_different_allele_is_a_false_positive_and_a_miss():
    truth = {("c", 10, "A", "G"): "SNP"}
    calls = {("c", 10, "A", "T"): "SNP"}

    stats = compare_variants(truth, calls)

    assert stats["SNP"]["tp"] == 0
    assert stats["SNP"]["fn"] == 1
    assert stats["SNP"]["fp"] == 1


VCF_WITH_CLASS = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "c\t10\t.\tA\tG\t.\tPASS\tSAMPLE=s;CLASS=SNP;DP=9\n"
    "c\t40\t.\tAT\tA\t.\tPASS\tCLASS=SMALL_INDEL\n"
)


def test_load_vcf_uses_the_class_in_info_and_skips_header_lines(tmp_path):
    path = tmp_path / "x.vcf"
    path.write_text(VCF_WITH_CLASS)

    assert load_vcf(path) == {
        ("c", 10, "A", "G"): "SNP",
        ("c", 40, "AT", "A"): "SMALL_INDEL",
    }


def test_load_vcf_derives_the_class_from_allele_lengths_when_info_has_none(tmp_path):
    big = "C" + "A" * 49  # deletes 49 bases
    bigger = "C" + "A" * 50  # deletes 50 bases
    text = (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "c\t10\t.\tA\tG\t.\t.\tDP=9\n"
        f"c\t20\t.\t{big}\tC\t.\t.\t.\n"
        f"c\t30\t.\tC\tC{'T' * 50}\t.\t.\t.\n"
        f"c\t40\t.\t{bigger}\tC\t.\t.\t.\n"
    )
    path = tmp_path / "bcftools.vcf"
    path.write_text(text)

    classes = {k[1]: v for k, v in load_vcf(path).items()}

    assert classes == {10: "SNP", 20: "SMALL_INDEL", 30: "LARGE_INDEL", 40: "LARGE_INDEL"}


def test_load_vcf_reads_gzip_files(tmp_path):
    path = tmp_path / "x.vcf.gz"
    path.write_bytes(gzip.compress(VCF_WITH_CLASS.encode()))

    assert len(load_vcf(path)) == 2


def test_main_prints_a_tab_separated_table_with_one_row_per_class_and_a_total(tmp_path, capsys):
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    truth = tmp_path / "truth.vcf"
    truth.write_text(header + "c\t10\t.\tA\tG\t.\t.\tCLASS=SNP\nc\t20\t.\tA\tC\t.\t.\tCLASS=SNP\n")
    calls = tmp_path / "calls.vcf"
    calls.write_text(header + "c\t10\t.\tA\tG\t.\t.\t.\nc\t99\t.\tA\tG\t.\t.\t.\n")

    main(["--truth", str(truth), "--calls", str(calls)])

    rows = [line.split("\t") for line in capsys.readouterr().out.splitlines()]
    assert rows[0] == ["class", "truth", "tp", "fn", "calls", "fp", "sensitivity", "precision"]
    by_class = {r[0]: r for r in rows[1:]}
    assert by_class["SNP"] == ["SNP", "2", "1", "1", "2", "1", "0.5000", "0.5000"]
    assert by_class["SMALL_INDEL"] == ["SMALL_INDEL", "0", "0", "0", "0", "0", "NA", "NA"]
    assert by_class["ALL"] == ["ALL", "2", "1", "1", "2", "1", "0.5000", "0.5000"]


def _del(pos, size):
    return ("c", pos, "C" + "A" * size, "C")


def _ins(pos, size):
    return ("c", pos, "C", "C" + "T" * size)


def _large(truth_keys, call_keys, **kw):
    truth = {k: "LARGE_INDEL" for k in truth_keys}
    calls = {k: "LARGE_INDEL" for k in call_keys}
    return compare_variants(truth, calls, **kw)["LARGE_INDEL"]


def test_a_large_call_that_is_a_few_bases_off_in_position_and_size_matches_the_truth():
    stats = _large([_del(5000, 100)], [_del(5012, 95)])

    assert (stats["tp"], stats["fn"], stats["fp"]) == (1, 0, 0)


def test_a_large_call_too_far_away_or_the_wrong_size_or_the_wrong_type_does_not_match():
    truth = [_del(5000, 100)]

    assert _large(truth, [_del(5080, 100)])["tp"] == 0   # 80 bp away, limit is 50
    assert _large(truth, [_del(5000, 60)])["tp"] == 0    # 60/100 is below 0.8
    assert _large(truth, [_ins(5000, 100)])["tp"] == 0   # insertion, not deletion


def test_two_large_calls_near_one_truth_variant_give_one_true_positive_and_one_false_positive():
    stats = _large([_del(5000, 100)], [_del(5003, 100), _del(5010, 98)])

    assert (stats["tp"], stats["fn"], stats["fp"]) == (1, 0, 1)


def test_small_indels_and_snps_still_need_an_exact_match():
    truth = {("c", 100, "AT", "A"): "SMALL_INDEL", ("c", 200, "A", "G"): "SNP"}
    calls = {("c", 101, "AT", "A"): "SMALL_INDEL", ("c", 201, "A", "G"): "SNP"}

    stats = compare_variants(truth, calls)

    assert stats["SMALL_INDEL"]["tp"] == 0
    assert stats["SNP"]["tp"] == 0


def test_a_symbolic_insertion_call_matches_an_explicit_truth_insertion_by_its_svlen():
    truth = {_ins(5000, 100): "LARGE_INDEL"}
    sym = ("c", 5003, "C", "<INS>")
    calls = {sym: "LARGE_INDEL"}

    stats = compare_variants(truth, calls, call_svlen={sym: 97})["LARGE_INDEL"]

    assert (stats["tp"], stats["fn"], stats["fp"]) == (1, 0, 0)


def test_load_svlens_reads_svlen_from_info_for_symbolic_records_only(tmp_path):
    path = tmp_path / "x.vcf"
    path.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "c\t10\t.\tA\t<INS>\t.\tPASS\tCLASS=LARGE_INDEL;SVTYPE=INS;SVLEN=120;IMPRECISE\n"
        "c\t20\t.\tA\tG\t.\tPASS\tCLASS=SNP\n"
    )

    assert load_svlens(path) == {("c", 10, "A", "<INS>"): 120}


def test_main_matches_a_symbolic_large_insertion_call_against_an_explicit_truth_insertion(tmp_path, capsys):
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    truth = tmp_path / "truth.vcf"
    truth.write_text(header + f"c\t5000\t.\tC\tC{'T' * 100}\t.\t.\tCLASS=LARGE_INDEL\n")
    calls = tmp_path / "calls.vcf"
    calls.write_text(header + "c\t5004\t.\tC\t<INS>\t.\tPASS\tCLASS=LARGE_INDEL;SVTYPE=INS;SVLEN=96;IMPRECISE\n")

    main(["--truth", str(truth), "--calls", str(calls)])

    rows = {l.split("\t")[0]: l.split("\t") for l in capsys.readouterr().out.splitlines()}
    assert rows["LARGE_INDEL"] == ["LARGE_INDEL", "1", "1", "0", "1", "0", "1.0000", "1.0000"]

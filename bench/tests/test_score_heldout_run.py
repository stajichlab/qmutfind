from score_heldout_run import main, score_rows

HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _vcf(path, records):
    path.write_text(HEADER + "".join(f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\t{i}\n" for c, p, r, a, i in records))


def test_the_rows_count_truth_found_calls_and_false_calls_per_class(tmp_path):
    truth = tmp_path / "truth.vcf"
    calls = tmp_path / "calls.vcf"
    _vcf(truth, [("c1", 10, "A", "G", "CLASS=SNP"), ("c1", 20, "C", "T", "CLASS=SNP"), ("c1", 30, "AT", "A", "CLASS=SMALL_INDEL")])
    _vcf(calls, [("c1", 10, "A", "G", "AF=1"), ("c1", 99, "A", "G", "AF=1"), ("c1", 30, "AT", "A", "AF=1")])

    rows = {r["klass"]: r for r in score_rows(truth, calls)}

    assert (rows["SNP"]["truth"], rows["SNP"]["tp"], rows["SNP"]["calls"], rows["SNP"]["fp"]) == (2, 1, 2, 1)
    assert (rows["SMALL_INDEL"]["truth"], rows["SMALL_INDEL"]["tp"], rows["SMALL_INDEL"]["fp"]) == (1, 1, 0)
    assert rows["LARGE_INDEL"]["truth"] == 0


def test_only_truth_and_calls_on_the_listed_contigs_are_scored(tmp_path):
    truth = tmp_path / "truth.vcf"
    calls = tmp_path / "calls.vcf"
    _vcf(truth, [("c1", 10, "A", "G", "CLASS=SNP"), ("c2", 10, "A", "G", "CLASS=SNP")])
    _vcf(calls, [("c1", 10, "A", "G", "AF=1"), ("c2", 10, "A", "G", "AF=1"), ("c3", 5, "A", "G", "AF=1")])

    rows = {r["klass"]: r for r in score_rows(truth, calls, contigs={"c1", "c2"})}
    only_c1 = {r["klass"]: r for r in score_rows(truth, calls, contigs={"c1"})}

    assert (rows["SNP"]["truth"], rows["SNP"]["calls"]) == (2, 2)   # the call on c3 is outside the scored contigs
    assert (only_c1["SNP"]["truth"], only_c1["SNP"]["calls"]) == (1, 1)


def test_main_appends_rows_to_a_tsv_and_writes_the_header_once(tmp_path):
    truth = tmp_path / "truth.vcf"
    calls = tmp_path / "calls.vcf"
    _vcf(truth, [("c1", 10, "A", "G", "CLASS=SNP")])
    _vcf(calls, [("c1", 10, "A", "G", "AF=1")])
    out = tmp_path / "results.tsv"
    args = ["--truth", str(truth), "--calls", str(calls), "--out", str(out), "--seed", "3", "--mode", "standard", "--tech", "art", "--coverage", "30"]

    main(args)
    main(args)

    lines = out.read_text().splitlines()
    assert lines[0] == "seed\tmode\ttech\tcoverage\tklass\ttruth\ttp\tcalls\tfp"
    assert len(lines) == 1 + 2 * 3  # two runs, three classes each
    assert lines[1].split("\t")[:5] == ["3", "standard", "art", "30", "SNP"]

import gzip
import random
import io
from collections import Counter

from make_truth import (
    Variant,
    apply_variants,
    design_variants,
    main,
    read_fasta,
    write_fasta,
    write_truth_vcf,
)


def test_design_returns_the_requested_number_of_variants_in_each_class(reference):
    variants = design_variants(reference, n_snp=20, n_small=10, n_large=5, seed=1)

    counts = Counter(v.klass for v in variants)

    assert counts == {"SNP": 20, "SMALL_INDEL": 10, "LARGE_INDEL": 5}


def test_every_variant_ref_allele_matches_the_reference_at_its_position(reference):
    variants = design_variants(reference, n_snp=20, n_small=10, n_large=5, seed=1)

    for v in variants:
        start = v.pos - 1
        assert reference[v.contig][start : start + len(v.ref)] == v.ref, v


def test_variants_on_a_contig_are_at_least_200_bp_apart(reference):
    # Many variants on small contigs make crowding likely, so this test would catch a missing rule.
    variants = design_variants(reference, n_snp=100, n_small=40, n_large=15, seed=3)

    for contig in reference:
        spans = sorted(
            (v.pos, v.pos + len(v.ref) - 1) for v in variants if v.contig == contig
        )
        for (_, end_a), (start_b, _) in zip(spans, spans[1:]):
            assert start_b - end_a >= 200, (contig, end_a, start_b)


FASTA_TEXT = ">chrA some description\nACGTac\ngtNN\n>chrB\nTTTT\n"


def test_read_fasta_gives_the_same_uppercase_result_for_plain_and_gzip_input(tmp_path):
    plain = tmp_path / "ref.fa"
    plain.write_text(FASTA_TEXT)
    zipped = tmp_path / "ref.fa.gz"
    zipped.write_bytes(gzip.compress(FASTA_TEXT.encode()))

    expected = {"chrA": "ACGTACGTNN", "chrB": "TTTT"}

    assert read_fasta(plain) == expected
    assert read_fasta(zipped) == expected


def test_read_fasta_reads_zstd_input_and_fails_loudly_on_a_corrupt_file(tmp_path):
    import subprocess

    import pytest

    plain = tmp_path / "ref.fa"
    plain.write_text(FASTA_TEXT)
    subprocess.run(["zstd", "-q", str(plain), "-o", str(tmp_path / "ref.fa.zst")], check=True)
    assert read_fasta(tmp_path / "ref.fa.zst") == {"chrA": "ACGTACGTNN", "chrB": "TTTT"}

    bad = tmp_path / "bad.fa.zst"
    bad.write_bytes(b"this is not zstd data")
    with pytest.raises(RuntimeError):
        read_fasta(bad)


def test_write_fasta_wraps_at_80_columns_and_round_trips(reference, tmp_path):
    path = tmp_path / "out.fa"
    with open(path, "w") as fh:
        write_fasta(reference, fh)

    lines = path.read_text().splitlines()
    assert all(len(l) <= 80 for l in lines)
    assert lines[0] == ">chrA"
    assert read_fasta(path) == reference


def _run_main(reference, tmp_path, prefix, seed):
    ref_path = tmp_path / "ref.fa"
    with open(ref_path, "w") as fh:
        write_fasta(reference, fh)
    out_prefix = tmp_path / prefix
    main([
        "--reference", str(ref_path), "--out-prefix", str(out_prefix), "--seed", str(seed),
        "--snps", "20", "--small-indels", "10", "--large-indels", "5",
    ])
    return out_prefix


def test_main_writes_a_truth_vcf_that_reproduces_the_mutated_fasta(reference, tmp_path):
    out_prefix = _run_main(reference, tmp_path, "run1", seed=1)

    variants = []
    for line in (tmp_path / "run1.truth.vcf").read_text().splitlines():
        if line.startswith("#"):
            continue
        chrom, pos, _id, ref, alt, _q, _f, info = line.split("\t")
        variants.append(Variant(chrom, int(pos), ref, alt, info.removeprefix("CLASS=")))

    assert len(variants) == 35
    assert apply_variants(reference, variants) == read_fasta(tmp_path / "run1.mutated.fa")


def test_main_with_the_same_seed_writes_identical_files(reference, tmp_path):
    _run_main(reference, tmp_path, "a", seed=9)
    _run_main(reference, tmp_path, "b", seed=9)

    assert (tmp_path / "a.truth.vcf").read_bytes() == (tmp_path / "b.truth.vcf").read_bytes()
    assert (tmp_path / "a.mutated.fa").read_bytes() == (tmp_path / "b.mutated.fa").read_bytes()


def test_truth_vcf_has_a_header_and_one_sorted_record_per_variant(reference):
    variants = design_variants(reference, n_snp=20, n_small=10, n_large=5, seed=1)
    out = io.StringIO()

    write_truth_vcf(variants, reference, out)

    lines = out.getvalue().splitlines()
    assert lines[0] == "##fileformat=VCFv4.2"
    assert "##contig=<ID=chrA,length=50000>" in lines
    assert "##contig=<ID=chrB,length=30000>" in lines
    assert any(l.startswith("##INFO=<ID=CLASS,") for l in lines)
    columns = [l for l in lines if l.startswith("#CHROM")]
    assert columns == ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    records = [l.split("\t") for l in lines if not l.startswith("#")]
    assert len(records) == 35
    order = {"chrA": 0, "chrB": 1}
    keys = [(order[r[0]], int(r[1])) for r in records]
    assert keys == sorted(keys)
    by_site = {(v.contig, v.pos): v for v in variants}
    for chrom, pos, _id, ref, alt, _qual, _filter, info in records:
        v = by_site[(chrom, int(pos))]
        assert (ref, alt, info) == (v.ref, v.alt, f"CLASS={v.klass}")


def test_no_variant_touches_a_base_outside_acgt(reference):
    gapped = dict(reference)
    seq = gapped["chrA"]
    gapped["chrA"] = seq[:10_000] + "N" * 20_000 + seq[30_000:]

    variants = design_variants(gapped, n_snp=60, n_small=20, n_large=5, seed=5)

    for v in variants:
        start = v.pos - 1
        footprint = gapped[v.contig][start : start + len(v.ref)]
        assert set(footprint) <= set("ACGT"), v


def test_applying_variants_changes_each_contig_length_by_inserted_minus_deleted(reference):
    variants = design_variants(reference, n_snp=20, n_small=10, n_large=5, seed=1)

    mutated = apply_variants(reference, variants)

    for contig, seq in reference.items():
        delta = sum(
            len(v.alt) - len(v.ref) for v in variants if v.contig == contig
        )
        assert len(mutated[contig]) == len(seq) + delta


def _positions_by_contig(variants):
    out = {}
    for v in variants:
        out.setdefault(v.contig, []).append((v.pos, v.pos + len(v.ref) - 1))
    return {c: sorted(spans) for c, spans in out.items()}


def test_min_spacing_is_a_setting_and_a_smaller_value_allows_closer_variants():
    ref = {"c": "".join(random.Random(3).choice("ACGT") for _ in range(8_000))}

    variants = design_variants(ref, n_snp=40, n_small=0, n_large=0, seed=1, min_spacing=30)

    spans = _positions_by_contig(variants)["c"]
    gaps = [b[0] - a[1] for a, b in zip(spans, spans[1:])]
    assert min(gaps) >= 30
    assert min(gaps) < 200  # with 200 bp spacing the smallest gap could not be this small


def test_homopolymer_fraction_puts_that_share_of_small_indels_in_runs_of_four_or_more_bases():
    rng = random.Random(5)
    ref = {"c": "".join(rng.choice("ACGT") for _ in range(120_000))}

    variants = design_variants(ref, n_snp=0, n_small=40, n_large=0, seed=2, homopolymer_fraction=1.0)

    seq = ref["c"]
    assert len(variants) == 40
    for v in variants:
        assert v.klass == "SMALL_INDEL"
        assert abs(len(v.ref) - len(v.alt)) == 1  # one base
        # The changed base is the run: seq[pos:pos+4] (0-based, after the anchor) is four copies of it.
        changed = v.ref[1:] if len(v.ref) > len(v.alt) else v.alt[1:]
        assert seq[v.pos : v.pos + 4] == changed * 4, v


def test_a_homopolymer_fraction_below_one_leaves_the_rest_of_the_small_indels_ordinary():
    rng = random.Random(6)
    ref = {"c": "".join(rng.choice("ACGT") for _ in range(120_000))}

    variants = design_variants(ref, n_snp=0, n_small=60, n_large=0, seed=3, homopolymer_fraction=0.5)

    in_run = 0
    for v in variants:
        changed = v.ref[1:] if len(v.ref) > len(v.alt) else v.alt[1:]
        if len(changed) == 1 and ref["c"][v.pos : v.pos + 4] == changed * 4:
            in_run += 1
    assert 30 <= in_run <= 45  # 30 are placed in runs by design; a few ordinary indels land in one by chance


def test_cluster_fraction_places_variants_close_to_an_earlier_one_without_overlapping_it():
    rng = random.Random(7)
    ref = {"c": "".join(rng.choice("ACGT") for _ in range(200_000))}

    variants = design_variants(ref, n_snp=80, n_small=20, n_large=0, seed=4, cluster_fraction=0.5)

    spans = _positions_by_contig(variants)["c"]
    close = 0
    for i, (s, e) in enumerate(spans):
        assert i == 0 or s > spans[i - 1][1], "variants overlap"
        near_before = i > 0 and s - spans[i - 1][1] <= 40
        near_after = i + 1 < len(spans) and spans[i + 1][0] - e <= 40
        close += near_before or near_after
    assert close >= 40  # half of the 100 variants were placed within 40 bp of another, so at least this many are close

from make_genome_list import check_entry, make_entry, parse_summary

HEADER = (
    "#assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\ttaxid\tspecies_taxid\torganism_name\t"
    "infraspecific_name\tisolate\tversion_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\tasm_name\t"
    "asm_submitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\texcluded_from_refseq\trelation_to_type_material\t"
    "asm_not_live_date\tassembly_type\tgroup\tgenome_size\tgenome_size_ungapped\tgc_percent\treplicon_count\t"
    "scaffold_count\tcontig_count\tannotation_provider\tannotation_name\tannotation_date\ttotal_gene_count\t"
    "protein_coding_gene_count\tnon_coding_gene_count\tpubmed_id"
)


def row(accession="GCF_000002655.1", organism="Aspergillus fumigatus Af293", strain="strain=Af293", status="latest",
        level="Chromosome", asm="ASM265v1", ftp="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1",
        size="29384958", provider="VEuPathDB", genes="10086", excluded=""):
    f = [accession, "PRJ", "SAM", "na", "reference genome", "1", "2", organism, strain, "na", status, level, "Major",
         "Full", "2007/04/03", asm, "sub", "GCA_x", "different", ftp, excluded, "na", "na", "haploid", "fungi", size,
         size, "50.0", "8", "8", "19", provider, "Name", "2026-10-19", genes, "9823", "262", "1"]
    return "\t".join(f)


def table(*rows):
    return "##comment\n" + HEADER + "\n" + "\n".join(rows) + "\n"


def test_the_table_is_read_into_rows_keyed_by_accession_and_the_comment_line_is_skipped():
    rows = parse_summary(table(row(), row(accession="GCF_000000001.1", organism="Other fungus")))

    assert set(rows) == {"GCF_000002655.1", "GCF_000000001.1"}
    assert rows["GCF_000002655.1"]["organism_name"] == "Aspergillus fumigatus Af293"


def test_an_entry_has_https_urls_for_the_fasta_and_the_gff_built_from_the_ftp_path():
    entry = make_entry(parse_summary(table(row()))["GCF_000002655.1"])

    base = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1"
    assert entry["fasta"]["url"] == base + "_genomic.fna.gz"
    assert entry["gff"]["url"] == base + "_genomic.gff.gz"
    assert entry["fasta"]["name"] == "GCF_000002655.1_ASM265v1_genomic.fna.gz"
    assert entry["accession"] == "GCF_000002655.1" and entry["genome_size"] == 29384958 and entry["genes"] == 10086


def test_an_ftp_address_is_turned_into_https():
    r = parse_summary(table(row(ftp="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1")))
    entry = make_entry(r["GCF_000002655.1"])

    assert entry["fasta"]["url"].startswith("https://ftp.ncbi.nlm.nih.gov/")


def test_the_label_adds_the_strain_only_when_the_organism_name_does_not_have_it():
    r = parse_summary(table(row(), row(accession="GCF_000000002.1", organism="Nakaseomyces glabratus", strain="strain=ATCC 2001"),
                            row(accession="GCF_000000003.1", organism="Candida sp", strain="na")))

    assert make_entry(r["GCF_000002655.1"])["label"] == "Aspergillus fumigatus Af293"
    assert make_entry(r["GCF_000000002.1"])["label"] == "Nakaseomyces glabratus ATCC 2001"
    assert make_entry(r["GCF_000000003.1"])["label"] == "Candida sp"


def test_a_genome_without_annotation_a_withdrawn_one_or_a_huge_one_is_rejected_with_the_reason():
    ok = make_entry(parse_summary(table(row()))["GCF_000002655.1"])
    assert check_entry(ok, max_size=100_000_000) == []

    no_annotation = make_entry(parse_summary(table(row(provider="na", genes="0")))["GCF_000002655.1"])
    assert any("annotation" in e for e in check_entry(no_annotation, max_size=100_000_000))

    old = make_entry(parse_summary(table(row(status="suppressed")))["GCF_000002655.1"])
    assert any("latest" in e for e in check_entry(old, max_size=100_000_000))

    huge = make_entry(parse_summary(table(row(size="150000000")))["GCF_000002655.1"])
    assert any("size" in e for e in check_entry(huge, max_size=100_000_000))


def test_a_trailing_slash_in_the_ftp_path_does_not_double_the_slash_in_the_urls():
    # The real table writes the path with a trailing slash.
    r = parse_summary(table(row(ftp="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/")))
    entry = make_entry(r["GCF_000002655.1"])

    assert "//_genomic" not in entry["fasta"]["url"] and "ASM265v1//" not in entry["gff"]["url"]
    assert entry["fasta"]["url"].endswith("/GCF_000002655.1_ASM265v1/GCF_000002655.1_ASM265v1_genomic.fna.gz")
    assert entry["fasta"]["name"] == "GCF_000002655.1_ASM265v1_genomic.fna.gz"


def test_selection_takes_every_chromosome_or_complete_genome_up_to_the_size_limit_plus_the_named_ones():
    from make_genome_list import select_accessions

    rows = parse_summary(table(
        row(accession="GCF_000000001.1", level="Chromosome", size="30000000"),
        row(accession="GCF_000000002.1", level="Complete Genome", size="20000000"),
        row(accession="GCF_000000003.1", level="Scaffold", size="25000000"),
        row(accession="GCF_000000004.1", level="Contig", size="25000000"),
        row(accession="GCF_000000005.1", level="Chromosome", size="90000000"),  # too large
        row(accession="GCF_000000006.1", level="Chromosome", status="suppressed"),
        row(accession="GCF_000000007.1", level="Chromosome", provider="na", genes="0"),  # no annotation
    ))

    chosen = select_accessions(rows, max_size=70_000_000, named=["GCF_000000003.1", "GCF_000000001.1"])

    assert chosen == ["GCF_000000001.1", "GCF_000000002.1", "GCF_000000003.1"]  # sorted, no duplicate, scaffold only by name


def test_a_genbank_entry_has_no_gff_when_ncbi_has_none_but_a_refseq_entry_without_one_is_a_problem(monkeypatch):
    import urllib.error

    import make_genome_list as m

    def fake_head(url):
        if url.endswith(".gff.gz"):
            raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)
        return 200, 1000, "*"

    monkeypatch.setattr(m, "head", fake_head)
    genbank = m.make_entry(parse_summary(table(row(accession="GCA_000000009.1", provider="na", genes="0")))["GCA_000000009.1"])
    refseq = m.make_entry(parse_summary(table(row()))["GCF_000002655.1"])

    assert genbank["source"] == "GenBank" and refseq["source"] == "RefSeq"
    assert m.add_file_sizes(genbank) == [] and genbank["gff"] is None and genbank["fasta"]["bytes"] == 1000
    assert any("gff" in p for p in m.add_file_sizes(refseq))


def test_a_busy_server_answer_503_or_a_timeout_is_retried_but_a_404_is_not(monkeypatch):
    import urllib.error

    import make_genome_list as m

    monkeypatch.setattr(m.time, "sleep", lambda s: None)
    calls = {"n": 0}

    def flaky(url):
        calls["n"] += 1
        if calls["n"] < 3:
            raise urllib.error.HTTPError(url, 503, "Service Unavailable", {}, None)
        return 200, 5, "*"

    monkeypatch.setattr(m, "head_once", flaky)
    assert m.head("https://x/y") == (200, 5, "*") and calls["n"] == 3

    def missing(url):
        calls["n"] += 1
        raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)

    calls["n"] = 0
    monkeypatch.setattr(m, "head_once", missing)
    try:
        m.head("https://x/y")
        assert False, "expected HTTPError"
    except urllib.error.HTTPError:
        pass
    assert calls["n"] == 1


def test_a_socket_timeout_is_retried(monkeypatch):
    import socket

    import make_genome_list as m

    monkeypatch.setattr(m.time, "sleep", lambda s: None)
    calls = {"n": 0}

    def slow(url):
        calls["n"] += 1
        if calls["n"] < 2:
            raise socket.timeout("The read operation timed out")
        return 200, 7, "*"

    monkeypatch.setattr(m, "head_once", slow)
    assert m.head("https://x/y") == (200, 7, "*")

from compare_annotation import (
    compare_annotations,
    parse_snpeff_ann,
    effect_matches,
    parse_bcsq,
    parse_csq_aa,
    parse_my_aa,
)


def test_parse_bcsq_reads_every_entry_and_splits_combined_consequences():
    info = "CLASS=SNP;BCSQ=frameshift&splice_region|G1|tx1|protein_coding|+|10A>10fs|100AC>A,@100"

    entries = parse_bcsq(info)

    assert len(entries) == 1  # the "@100" entry only points at another record and is dropped
    e = entries[0]
    assert e["terms"] == ["frameshift", "splice_region"]
    assert (e["gene"], e["transcript"], e["strand"]) == ("G1", "tx1", "+")


def test_a_record_without_bcsq_has_no_entries():
    assert parse_bcsq("DP=20;AF=1.0") == []


def test_amino_acid_changes_are_parsed_from_both_tools_into_the_same_form():
    assert parse_csq_aa("230V>230A") == ("V", 230, "A")
    assert parse_csq_aa("565G") == ("G", 565, "G")        # synonymous: one residue
    assert parse_csq_aa("12*>12Q") == ("*", 12, "Q")
    assert parse_csq_aa("463YGNF>463Y") is None            # an indel is not a single change
    assert parse_my_aa("A2D") == ("A", 2, "D")
    assert parse_my_aa("W9*") == ("W", 9, "*")
    assert parse_my_aa("*17Q") == ("*", 17, "Q")
    assert parse_my_aa(None) is None


def test_effect_matching_maps_my_labels_to_sequence_ontology_terms():
    assert effect_matches("missense", ["missense"])
    assert effect_matches("nonsense", ["stop_gained"])
    assert effect_matches("synonymous", ["stop_retained"])
    assert effect_matches("frameshift", ["frameshift", "splice_region"])
    assert effect_matches("splice_site", ["splice_donor"]) and effect_matches("splice_site", ["splice_acceptor"])
    assert effect_matches("intron_derived", ["intron"])
    assert not effect_matches("missense", ["synonymous"])
    assert not effect_matches("frameshift", ["inframe_deletion"])


def test_an_entry_with_no_transcript_is_matched_to_my_transcript_by_gene():
    mine = {("c", 5, "A", "C"): {"effect": "intron_derived", "transcript": "tx1", "gene": "G1", "aa": None}}
    csq = {("c", 5, "A", "C"): [{"terms": ["intron"], "gene": "G1", "transcript": "", "aa": ""}]}

    r = compare_annotations(mine, csq)

    assert (r["agree"], r["csq_other_transcript"]) == (1, 0)


def test_csq_calls_the_utr_an_intron_because_it_has_no_utr_features_and_that_is_not_a_disagreement():
    mine = {("c", 5, "A", "C"): {"effect": "5UTR_derived", "transcript": "tx1", "gene": "G1", "aa": None}}
    csq = {("c", 5, "A", "C"): [{"terms": ["intron"], "gene": "G1", "transcript": "", "aa": ""}]}

    r = compare_annotations(mine, csq)

    assert (r["disagree"], r["utr_vs_csq_intron"]) == (0, 1)


def test_a_gene_written_as_its_name_by_csq_is_matched_to_my_locus_tag_with_an_alias_map():
    mine = {("c", 5, "A", "C"): {"effect": "intron_derived", "transcript": "tx1", "gene": "AFUA_1G03150", "aa": None}}
    csq = {("c", 5, "A", "C"): [{"terms": ["intron"], "gene": "erg24", "transcript": "", "aa": ""}]}

    without = compare_annotations(mine, csq)
    with_aliases = compare_annotations(mine, csq, gene_aliases={"AFUA_1G03150": {"erg24"}})

    assert without["csq_other_transcript"] == 1
    assert (with_aliases["agree"], with_aliases["csq_other_transcript"]) == (1, 0)


MINE = {
    ("c", 10, "A", "G"): {"effect": "missense", "transcript": "tx1", "aa": "A2D"},
    ("c", 20, "C", "T"): {"effect": "synonymous", "transcript": "tx1", "aa": "E3E"},
    ("c", 30, "G", "A"): {"effect": "nonsense", "transcript": "tx1", "aa": "W9*"},
    ("c", 40, "T", "A"): {"effect": "intergenic", "transcript": None, "aa": None},
    ("c", 50, "A", "C"): {"effect": "5UTR_derived", "transcript": "tx1", "aa": None},
    ("c", 60, "A", "C"): {"effect": "missense", "transcript": "tx1", "aa": "K5R"},
}
CSQ = {
    ("c", 10, "A", "G"): [{"terms": ["missense"], "transcript": "tx1", "aa": "2A>2D"}],
    ("c", 20, "C", "T"): [{"terms": ["missense"], "transcript": "tx1", "aa": "3E>3K"}],   # disagrees
    ("c", 30, "G", "A"): [{"terms": ["stop_gained"], "transcript": "tx2", "aa": "9W>9*"}],  # other transcript only
    ("c", 40, "T", "A"): [],                                                               # both silent
    ("c", 50, "A", "C"): [],                                                               # csq has no UTR features
    ("c", 60, "A", "C"): [{"terms": ["missense"], "transcript": "tx1", "aa": "5K>5Q"}],    # same effect, other amino acid
}


def test_the_comparison_counts_agreement_disagreement_and_silence_on_my_transcript():
    r = compare_annotations(MINE, CSQ)

    assert r["total"] == 6
    assert r["agree"] == 2                    # 10 and 60 have a matching term on tx1
    assert r["disagree"] == 1                 # 20: synonymous against missense
    assert r["csq_other_transcript"] == 1     # 30: csq reports only tx2
    assert r["both_silent"] == 1              # 40: intergenic in mine, nothing in csq
    assert r["csq_silent_mine_genic"] == 1    # 50: mine says UTR, csq says nothing


def test_the_amino_acid_comparison_covers_coding_snps_that_both_tools_report_on_my_transcript():
    r = compare_annotations(MINE, CSQ)

    assert r["aa_compared"] == 3              # 10, 20 and 60
    assert r["aa_agree"] == 1                 # only 10 has the same residues and position
    assert r["aa_disagree"] == 2


def test_the_examples_list_the_disagreements_with_both_tools_calls():
    r = compare_annotations(MINE, CSQ)

    kinds = {e["kind"] for e in r["examples"]}
    assert "disagree" in kinds and "aa_disagree" in kinds
    first = next(e for e in r["examples"] if e["kind"] == "disagree")
    assert first["key"] == ("c", 20, "C", "T") and first["mine"] == "synonymous"


ANN_MISSENSE = "G|missense_variant|MODERATE|g1|GENE_x|transcript|rna-XM_1.1|protein_coding|1/1|c.689T>C|p.Val230Ala|689/1197|689/1197|230/398||"
ANN_UPSTREAM = "G|upstream_gene_variant|MODIFIER|g2|gene-g2|transcript|XM_2.1|protein_coding||c.-4462A>G|||||4462|"


def test_snpeff_entries_use_my_term_names_the_transcript_id_and_the_same_amino_acid_form():
    entries = parse_snpeff_ann(f"CLASS=SNP;ANN={ANN_MISSENSE},{ANN_UPSTREAM}")

    assert len(entries) == 1  # the upstream entry says the variant is outside the gene and is dropped
    e = entries[0]
    assert (e["terms"], e["transcript"], e["aa"]) == (["missense"], "rna-XM_1.1", "230V>230A")


def test_snpeff_terms_are_normalised_to_the_names_the_comparison_uses():
    def terms(annotation):
        return parse_snpeff_ann(f"ANN=T|{annotation}|HIGH|g|G|transcript|rna-X|protein_coding|1/1|c.1A>T||||||")[0]["terms"]

    assert terms("stop_gained") == ["stop_gained"]
    assert terms("frameshift_variant") == ["frameshift"]
    assert terms("splice_donor_variant&intron_variant") == ["splice_donor", "intron"]
    assert terms("5_prime_UTR_variant") == ["5_prime_utr"]
    assert terms("3_prime_UTR_variant") == ["3_prime_utr"]
    assert terms("disruptive_inframe_insertion") == ["inframe_insertion"]
    assert terms("conservative_inframe_deletion") == ["inframe_deletion"]
    assert terms("synonymous_variant") == ["synonymous"]
    assert terms("non_coding_transcript_exon_variant") == ["non_coding"]


def test_snpeff_protein_changes_are_converted_from_three_letter_hgvs():
    def aa(hgvs):
        return parse_snpeff_ann(f"ANN=T|missense_variant|MODERATE|g|G|transcript|rna-X|protein_coding|1/1|c.1A>T|{hgvs}||||||")[0]["aa"]

    assert aa("p.Val230Ala") == "230V>230A"
    assert aa("p.Trp9*") == "9W>9*"
    assert aa("p.Glu3Glu") == "3E"          # synonymous: one residue, as bcftools csq writes it
    assert aa("p.Glu11fs") == ""             # a frameshift has no single change
    assert aa("") == ""


def test_a_record_without_ann_has_no_snpeff_entries():
    assert parse_snpeff_ann("DP=20") == []

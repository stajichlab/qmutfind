//! Effect of a variant on the canonical transcript. The fixtures are small and hand-checked.
//!
//! Contig c1 (200 bp), a plus-strand gene with two exons and one intron:
//!   21-23 5' UTR, 24-50 CDS (codons 1-9), 51-80 intron (GT...AG), 81-104 CDS (codons 10-17, last is TAA), 105-110 3' UTR
//!   codons: ATG GCT GAA CGT AAA TTA GAT CCC TGG | GGC AAA CTG GAC CAT GAA TGC TAA  (M A E R K L D P W G K L D H E C *)
//! Contig c2 (40 bp), a minus-strand gene, one exon 11-22, CDS 11-22, coding strand ATG AAA CCC TAA (M K P *).
//! Contig c3 (60 bp), a minus-strand gene, exon 11-40, CDS 20-31: forward 32-40 is the 5' UTR, 11-19 the 3' UTR.
//! Contig c4 (60 bp), a plus-strand non-coding RNA, exons 11-20 and 31-40, intron 21-30.
//! Contig c5 (80 bp), two overlapping plus-strand genes. P: exon 11-60, CDS 50-58 (so 11-49 is its 5' UTR).
//! Q: exon 25-45, CDS 28-42 (ATG GCT GAA CGT AAA). Position 35 is in P's UTR and in the CDS of Q.

use qmutfind_core::{annotate_text, annotate_vcf, classify_effect, Annotation, Genome};

fn c1() -> Vec<u8> {
    let mut s = String::new();
    s.push_str(&"T".repeat(20)); // 1-20
    s.push_str("CAC"); // 21-23 5' UTR
    s.push_str("ATGGCTGAACGTAAATTAGATCCCTGG"); // 24-50
    s.push_str(&format!("GT{}AG", "A".repeat(26))); // 51-80 intron
    s.push_str("GGCAAACTGGACCATGAATGCTAA"); // 81-104
    s.push_str("CCCCCC"); // 105-110 3' UTR
    s.push_str(&"T".repeat(90)); // 111-200
    assert_eq!(s.len(), 200);
    s.into_bytes()
}

fn c2() -> Vec<u8> {
    let s = format!("{}TTAGGGTTTCAT{}", "T".repeat(10), "T".repeat(18));
    assert_eq!(s.len(), 40);
    s.into_bytes()
}

fn c5() -> Vec<u8> {
    // 1-27 filler, 28-42 CDS of Q, 43-49 filler, 50-58 CDS of P (ATGAAAAAA), 59-80 filler
    let s = format!("{}ATGGCTGAACGTAAA{}ATGAAAAAA{}", "T".repeat(27), "T".repeat(7), "T".repeat(22));
    assert_eq!(s.len(), 80);
    s.into_bytes()
}

const GFF: &str = "##gff-version 3
c1\tsrc\tgene\t21\t110\t.\t+\t.\tID=gene1;locus_tag=G1
c1\tsrc\tmRNA\t21\t110\t.\t+\t.\tID=tx1;Parent=gene1
c1\tsrc\texon\t21\t50\t.\t+\t.\tParent=tx1
c1\tsrc\texon\t81\t110\t.\t+\t.\tParent=tx1
c1\tsrc\tCDS\t24\t50\t.\t+\t0\tParent=tx1
c1\tsrc\tCDS\t81\t104\t.\t+\t0\tParent=tx1
c2\tsrc\tgene\t11\t22\t.\t-\t.\tID=gene2;locus_tag=G2
c2\tsrc\tmRNA\t11\t22\t.\t-\t.\tID=tx2;Parent=gene2
c2\tsrc\texon\t11\t22\t.\t-\t.\tParent=tx2
c2\tsrc\tCDS\t11\t22\t.\t-\t0\tParent=tx2
c3\tsrc\tgene\t11\t40\t.\t-\t.\tID=gene3;locus_tag=G3
c3\tsrc\tmRNA\t11\t40\t.\t-\t.\tID=tx3;Parent=gene3
c3\tsrc\texon\t11\t40\t.\t-\t.\tParent=tx3
c3\tsrc\tCDS\t20\t31\t.\t-\t0\tParent=tx3
c4\tsrc\tgene\t11\t40\t.\t+\t.\tID=gene4;locus_tag=G4
c4\tsrc\tncRNA\t11\t40\t.\t+\t.\tID=tx4;Parent=gene4
c4\tsrc\texon\t11\t20\t.\t+\t.\tParent=tx4
c4\tsrc\texon\t31\t40\t.\t+\t.\tParent=tx4
c5\tsrc\tgene\t11\t60\t.\t+\t.\tID=geneP;locus_tag=GP
c5\tsrc\tmRNA\t11\t60\t.\t+\t.\tID=txP;Parent=geneP
c5\tsrc\texon\t11\t60\t.\t+\t.\tParent=txP
c5\tsrc\tCDS\t50\t58\t.\t+\t0\tParent=txP
c5\tsrc\tgene\t25\t45\t.\t+\t.\tID=geneQ;locus_tag=GQ
c5\tsrc\tmRNA\t25\t45\t.\t+\t.\tID=txQ;Parent=geneQ
c5\tsrc\texon\t25\t45\t.\t+\t.\tParent=txQ
c5\tsrc\tCDS\t28\t42\t.\t+\t0\tParent=txQ
";

fn setup() -> (Annotation, Genome) {
    let ann = Annotation::from_gff(GFF.as_bytes()).unwrap();
    let genome = Genome::new(vec![
        ("c1".to_string(), c1()),
        ("c2".to_string(), c2()),
        ("c3".to_string(), "T".repeat(60).into_bytes()),
        ("c4".to_string(), "T".repeat(60).into_bytes()),
        ("c5".to_string(), c5()),
    ]);
    (ann, genome)
}

fn effect(contig: &str, pos: u64, reference: &str, alt: &str) -> (String, Option<String>) {
    let (ann, genome) = setup();
    let e = classify_effect(&ann, &genome, contig, pos, reference, alt);
    (e.label, e.aa)
}

#[test]
fn the_gff_gives_one_canonical_transcript_per_gene_with_the_locus_tag_and_transcript_id() {
    let (ann, genome) = setup();

    let e = classify_effect(&ann, &genome, "c1", 28, "C", "A");

    assert_eq!((e.gene.as_deref(), e.transcript.as_deref(), e.n_transcripts), (Some("G1"), Some("tx1"), 1));
}

#[test]
fn a_snp_that_keeps_the_amino_acid_is_synonymous() {
    // Codon 3 is GAA at 30-32. A to G at 32 gives GAG, still E.
    assert_eq!(effect("c1", 32, "A", "G"), ("synonymous".to_string(), Some("E3E".to_string())));
}

#[test]
fn a_snp_that_changes_the_amino_acid_is_missense_with_the_protein_change() {
    // Codon 2 is GCT at 27-29. C to A at 28 gives GAT, A to D.
    assert_eq!(effect("c1", 28, "C", "A"), ("missense".to_string(), Some("A2D".to_string())));
}

#[test]
fn a_snp_that_makes_a_stop_is_nonsense() {
    // Codon 9 is TGG at 48-50. G to A at 49 gives TAG.
    assert_eq!(effect("c1", 49, "G", "A"), ("nonsense".to_string(), Some("W9*".to_string())));
}

#[test]
fn a_snp_in_the_start_codon_is_start_lost() {
    // ATG at 24-26. T to C at 25 gives ACG.
    assert_eq!(effect("c1", 25, "T", "C"), ("start_lost".to_string(), Some("M1T".to_string())));
}

#[test]
fn a_snp_in_the_stop_codon_that_removes_the_stop_is_stop_lost_and_the_codon_number_counts_across_the_intron() {
    // TAA at 102-104 is codon 17. T to C at 102 gives CAA, Q.
    assert_eq!(effect("c1", 102, "T", "C"), ("stop_lost".to_string(), Some("*17Q".to_string())));
}

#[test]
fn a_snp_in_the_second_exon_uses_the_codon_phase_carried_over_from_the_first_exon() {
    // Codon 10 is GGC at 81-83 and codon 11 is AAA at 84-86. A to G at 85 gives AGA, K to R.
    assert_eq!(effect("c1", 85, "A", "G"), ("missense".to_string(), Some("K11R".to_string())));
    // The third base of the same codon, A to G at 86, gives AAG, still K.
    assert_eq!(effect("c1", 86, "A", "G"), ("synonymous".to_string(), Some("K11K".to_string())));
}

#[test]
fn on_the_minus_strand_the_codon_is_read_from_the_reverse_complement() {
    // Forward 11-22 is TTAGGGTTTCAT, coding ATGAAACCCTAA. Codon 2 AAA is forward 17-19 (TTT).
    // T to C at 18 reads as A to G on the coding strand: AGA, K to R.
    assert_eq!(effect("c2", 18, "T", "C"), ("missense".to_string(), Some("K2R".to_string())));
}


#[test]
fn exon_bases_outside_the_cds_are_utrs_named_by_strand() {
    assert_eq!(effect("c1", 22, "T", "A"), ("5UTR_derived".to_string(), None));
    assert_eq!(effect("c1", 107, "C", "A"), ("3UTR_derived".to_string(), None));
    // On the minus strand the 5' side is the high-coordinate side.
    assert_eq!(effect("c3", 35, "T", "A"), ("5UTR_derived".to_string(), None));
    assert_eq!(effect("c3", 15, "T", "A"), ("3UTR_derived".to_string(), None));
}

#[test]
fn a_utr_variant_still_names_its_gene_and_transcript() {
    let (ann, genome) = setup();

    let e = classify_effect(&ann, &genome, "c1", 22, "T", "A");

    assert_eq!((e.gene.as_deref(), e.transcript.as_deref()), (Some("G1"), Some("tx1")));
}

#[test]
fn introns_are_derived_from_the_exon_gaps_and_the_two_bases_at_each_end_are_splice_sites() {
    for pos in [51, 52, 79, 80] {
        assert_eq!(effect("c1", pos, "A", "C").0, "splice_site", "position {pos}");
    }
    for pos in [53, 60, 78] {
        assert_eq!(effect("c1", pos, "A", "C").0, "intron_derived", "position {pos}");
    }
}

#[test]
fn a_snp_outside_every_transcript_or_on_a_contig_without_annotation_is_intergenic() {
    assert_eq!(effect("c1", 150, "T", "A"), ("intergenic".to_string(), None));
    assert_eq!(effect("c9", 10, "T", "A"), ("intergenic".to_string(), None));
}

#[test]
fn a_transcript_without_a_cds_is_non_coding_in_its_exons_and_has_derived_introns() {
    assert_eq!(effect("c4", 15, "T", "A").0, "non_coding");
    assert_eq!(effect("c4", 25, "T", "A").0, "intron_derived");
    assert_eq!(effect("c4", 21, "T", "A").0, "splice_site");
}

// ---- indels (VCF form: one anchor base, then the changed bases) ----

#[test]
fn an_indel_in_the_cds_is_a_frameshift_unless_its_length_is_a_multiple_of_three() {
    // c1 30-33 is G A A C. Delete A (31) or AAC (31-33), insert T or TTT after 30.
    assert_eq!(effect("c1", 30, "GA", "G").0, "frameshift");
    assert_eq!(effect("c1", 30, "GAAC", "G").0, "inframe_deletion");
    assert_eq!(effect("c1", 30, "G", "GT").0, "frameshift");
    // GGG keeps the reading frame and makes no stop. (TTT would give GCT GTT TAA: a new stop, tested below.)
    assert_eq!(effect("c1", 30, "G", "GGGG").0, "inframe_insertion");
    // Four and five bases are frameshifts too, not only one and two.
    assert_eq!(effect("c1", 30, "GAACG", "G").0, "frameshift");
    assert_eq!(effect("c1", 30, "G", "GTTTTT").0, "frameshift");
}

#[test]
fn a_frameshift_on_the_minus_strand_is_found_the_same_way() {
    // c2 forward 11-22 is TTAGGGTTTCAT. Delete the G at 16, anchored at 15.
    assert_eq!(effect("c2", 15, "GG", "G").0, "frameshift");
    assert_eq!(effect("c2", 15, "GGGT", "G").0, "inframe_deletion");
}

#[test]
fn an_indel_that_touches_a_splice_site_is_a_splice_site_variant() {
    // Delete the donor GT at 51-52, anchored at 50.
    assert_eq!(effect("c1", 50, "GGT", "G").0, "splice_site");
    // An insertion between the last exon base (50) and the first intron base (51).
    assert_eq!(effect("c1", 50, "G", "GA").0, "splice_site");
    // Delete the last exon base and the first intron base, anchored at 49.
    assert_eq!(effect("c1", 49, "GGG", "G").0, "splice_site");
}

#[test]
fn an_indel_wholly_inside_an_intron_or_a_utr_takes_the_label_of_that_region() {
    assert_eq!(effect("c1", 60, "AA", "A").0, "intron_derived");
    assert_eq!(effect("c1", 21, "CA", "C").0, "5UTR_derived");
    assert_eq!(effect("c1", 106, "CC", "C").0, "3UTR_derived");
}

#[test]
fn an_indel_that_crosses_from_the_cds_into_a_utr_is_complex() {
    // Delete 103-106: the end of the CDS (the stop codon) and the start of the 3' UTR.
    assert_eq!(effect("c1", 102, "TAACC", "T").0, "complex");
}

#[test]
fn a_symbolic_large_insertion_in_a_cds_has_no_known_frame() {
    assert_eq!(effect("c1", 30, "G", "<INS>").0, "large_indel_in_cds");
    assert_eq!(effect("c1", 60, "A", "<INS>").0, "intron_derived");
}

#[test]
fn an_indel_outside_every_transcript_is_intergenic_and_one_inside_names_its_gene() {
    assert_eq!(effect("c1", 150, "TT", "T").0, "intergenic");
    let (ann, genome) = setup();
    let e = classify_effect(&ann, &genome, "c1", 30, "GA", "G");
    assert_eq!((e.gene.as_deref(), e.transcript.as_deref()), (Some("G1"), Some("tx1")));
}


// ---- annotating a VCF ----

const VCF_HEADER: &str = "##fileformat=VCFv4.2\n##contig=<ID=c1,length=200>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn annotate(vcf_body: &str, with_annotation: bool) -> Vec<String> {
    let (ann, genome) = setup();
    let input = format!("{VCF_HEADER}{vcf_body}");
    let mut out = Vec::new();
    annotate_vcf(input.as_bytes(), &mut out, if with_annotation { Some(&ann) } else { None }, &genome).unwrap();
    String::from_utf8(out).unwrap().lines().map(str::to_string).collect()
}

#[test]
fn the_header_gets_the_five_info_definitions_before_the_column_line() {
    let lines = annotate("", true);

    assert_eq!(lines[0], "##fileformat=VCFv4.2");
    assert_eq!(lines[1], "##contig=<ID=c1,length=200>");
    for id in ["EFFECT", "GENE", "TRANSCRIPT", "AA", "NTX"] {
        assert!(lines.iter().any(|l| l.starts_with(&format!("##INFO=<ID={id},"))), "missing {id}");
    }
    assert_eq!(lines.last().unwrap(), "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
}

#[test]
fn a_record_keeps_its_info_and_gains_the_effect_tags() {
    let lines = annotate("c1\t28\t.\tC\tA\t.\tPASS\tDP=20;AF=1.0000\n", true);

    assert_eq!(
        lines.last().unwrap(),
        "c1\t28\t.\tC\tA\t.\tPASS\tDP=20;AF=1.0000;EFFECT=missense;GENE=G1;TRANSCRIPT=tx1;AA=A2D;NTX=1"
    );
}

#[test]
fn an_empty_info_field_is_replaced_not_extended() {
    let lines = annotate("c1\t150\t.\tT\tA\t.\tPASS\t.\n", true);

    assert_eq!(lines.last().unwrap(), "c1\t150\t.\tT\tA\t.\tPASS\tEFFECT=intergenic");
}

#[test]
fn without_an_annotation_every_record_is_unannotated() {
    let lines = annotate("c1\t28\t.\tC\tA\t.\tPASS\tDP=20\n", false);

    assert_eq!(lines.last().unwrap(), "c1\t28\t.\tC\tA\t.\tPASS\tDP=20;EFFECT=unannotated");
}

#[test]
fn a_multi_allele_record_gets_one_effect_per_allele_in_the_same_order() {
    // Position 28 is C in codon 2 (GCT). C to A gives GAT (D, missense). C to T gives GTT (V, missense).
    let lines = annotate("c1\t28\t.\tC\tA,T\t.\tPASS\t.\n", true);

    assert_eq!(
        lines.last().unwrap(),
        "c1\t28\t.\tC\tA,T\t.\tPASS\tEFFECT=missense,missense;GENE=G1;TRANSCRIPT=tx1;AA=A2D,A2V;NTX=1"
    );
}

#[test]
fn a_record_with_too_few_fields_is_an_error() {
    let (ann, genome) = setup();
    let input = format!("{VCF_HEADER}c1\t28\t.\tC\n");

    assert!(annotate_vcf(input.as_bytes(), &mut Vec::new(), Some(&ann), &genome).is_err());
}


#[test]
fn annotate_text_gives_the_same_output_as_annotate_vcf_with_and_without_a_gff() {
    let fasta = format!(">c1\n{}\n>c2\n{}\n", String::from_utf8(c1()).unwrap(), String::from_utf8(c2()).unwrap());
    let body = "c1\t28\t.\tC\tA\t.\tPASS\tDP=20\nc1\t150\t.\tT\tA\t.\tPASS\t.\n";
    let vcf = format!("{VCF_HEADER}{body}");

    let with_gff = annotate_text(&vcf, fasta.as_bytes(), Some(GFF.as_bytes())).unwrap();
    let without = annotate_text(&vcf, fasta.as_bytes(), None).unwrap();

    assert_eq!(with_gff.lines().map(str::to_string).collect::<Vec<_>>(), annotate(body, true));
    assert_eq!(without.lines().map(str::to_string).collect::<Vec<_>>(), annotate(body, false));
}

#[test]
fn annotate_text_reports_a_bad_gff_or_an_empty_reference_as_an_error() {
    let fasta = format!(">c1\n{}\n", String::from_utf8(c1()).unwrap());
    let vcf = format!("{VCF_HEADER}c1\t28\t.\tC\tA\t.\tPASS\t.\n");

    assert!(annotate_text(&vcf, fasta.as_bytes(), Some(b"c1\tsrc\tgene\tnot-a-number\t5\t.\t+\t.\tID=g\n")).is_err());
    assert!(annotate_text(&vcf, b"", Some(GFF.as_bytes())).is_err());
}


// ---- found by comparing with bcftools csq on 3,700 variants ----

#[test]
fn a_variant_in_the_utr_of_one_gene_and_the_cds_of_another_gets_the_more_severe_effect() {
    // c5 position 35 is in P's 5' UTR (P starts first) and in Q's codon 3 (GAA, 34-36). A to G gives GGA: E to G.
    let (ann, genome) = setup();

    let e = classify_effect(&ann, &genome, "c5", 35, "A", "G");

    assert_eq!((e.label.as_str(), e.aa.as_deref(), e.gene.as_deref()), ("missense", Some("E3G"), Some("GQ")));
}

#[test]
fn an_indel_whose_anchor_is_outside_the_transcript_is_judged_by_the_bases_it_affects() {
    // c1 18-20 are intergenic T, 21-23 are the 5' UTR CAC, 24 is the start of the CDS.
    // Delete 21-23, anchored at 20: every deleted base is in the UTR.
    assert_eq!(effect("c1", 20, "TCAC", "T").0, "5UTR_derived");
    // Delete 19-24, anchored at 18: two deleted bases are outside, three are UTR and one is CDS.
    assert_eq!(effect("c1", 18, "TTTCACA", "T").0, "complex");
}

#[test]
fn an_in_frame_insertion_or_deletion_that_gains_a_stop_loses_the_stop_or_loses_the_start_is_not_just_in_frame() {
    // Insert TAA between codon 3 (30-32) and codon 4: a new stop. Insert GGG: no stop.
    assert_eq!(effect("c1", 32, "A", "ATAA").0, "nonsense");
    assert_eq!(effect("c1", 32, "A", "AGGG").0, "inframe_insertion");
    // Delete the stop codon TAA (102-104), anchored at 101.
    assert_eq!(effect("c1", 101, "CTAA", "C").0, "stop_lost");
    // Delete the start codon ATG (24-26), anchored at 23.
    assert_eq!(effect("c1", 23, "CATG", "C").0, "start_lost");
    // Delete codon 16 (TGC, 99-101), anchored at 98: only the protein is shorter.
    assert_eq!(effect("c1", 98, "ATGC", "A").0, "inframe_deletion");
}

#[test]
fn an_in_frame_indel_that_gains_a_stop_on_the_minus_strand_is_found_the_same_way() {
    // c2 forward 11-22 is TTAGGGTTTCAT, coding ATGAAACCCTAA (M K P *), codon 2 is forward 17-19.
    // Insert the coding bases TAA after codon 1, which is the forward bases TTA between 19 and 20.
    // Forward insertion TTA after position 19 reads as TAA on the coding strand: a new stop.
    assert_eq!(effect("c2", 19, "T", "TTTA").0, "nonsense");
}

#[test]
fn a_deletion_anchored_on_the_last_base_of_a_gene_but_deleting_only_intergenic_bases_is_intergenic() {
    // c1 105-110 is the 3' UTR and the transcript ends at 110. Bases 111 and 112 are intergenic.
    // Delete 111-112, anchored at 110 (the anchor is in the gene, the deleted bases are not).
    assert_eq!(effect("c1", 110, "CTT", "C").0, "intergenic");
    // Delete 110-111 (anchored at 109): one deleted base is in the UTR and one is not, so it is complex.
    assert_eq!(effect("c1", 109, "CCT", "C").0, "complex");
}

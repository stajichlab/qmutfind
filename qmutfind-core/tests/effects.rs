//! Effect of a variant on the canonical transcript. The fixtures are small and hand-checked.
//!
//! Contig c1 (200 bp), a plus-strand gene with two exons and one intron:
//!   21-23 5' UTR, 24-50 CDS (codons 1-9), 51-80 intron (GT...AG), 81-104 CDS (codons 10-17, last is TAA), 105-110 3' UTR
//!   codons: ATG GCT GAA CGT AAA TTA GAT CCC TGG | GGC AAA CTG GAC CAT GAA TGC TAA  (M A E R K L D P W G K L D H E C *)
//! Contig c2 (40 bp), a minus-strand gene, one exon 11-22, CDS 11-22, coding strand ATG AAA CCC TAA (M K P *).
//! Contig c3 (60 bp), a minus-strand gene, exon 11-40, CDS 20-31: forward 32-40 is the 5' UTR, 11-19 the 3' UTR.
//! Contig c4 (60 bp), a plus-strand non-coding RNA, exons 11-20 and 31-40, intron 21-30.

use qmutfind_core::{classify_effect, Annotation, Genome};

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
";

fn setup() -> (Annotation, Genome) {
    let ann = Annotation::from_gff(GFF.as_bytes()).unwrap();
    let genome = Genome::new(vec![
        ("c1".to_string(), c1()),
        ("c2".to_string(), c2()),
        ("c3".to_string(), "T".repeat(60).into_bytes()),
        ("c4".to_string(), "T".repeat(60).into_bytes()),
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
    assert_eq!(effect("c1", 30, "G", "GTTT").0, "inframe_insertion");
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

//! Effect of a variant on the canonical transcript. The fixtures are small and hand-checked.
//!
//! Contig c1 (200 bp), a plus-strand gene with two exons and one intron:
//!   21-23 5' UTR, 24-50 CDS (codons 1-9), 51-80 intron (GT...AG), 81-104 CDS (codons 10-17, last is TAA), 105-110 3' UTR
//!   codons: ATG GCT GAA CGT AAA TTA GAT CCC TGG | GGC AAA CTG GAC CAT GAA TGC TAA  (M A E R K L D P W G K L D H E C *)
//! Contig c2 (40 bp), a minus-strand gene, one exon 11-22, CDS 11-22, coding strand ATG AAA CCC TAA (M K P *).

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
";

fn setup() -> (Annotation, Genome) {
    let ann = Annotation::from_gff(GFF.as_bytes()).unwrap();
    let genome = Genome::new(vec![("c1".to_string(), c1()), ("c2".to_string(), c2())]);
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

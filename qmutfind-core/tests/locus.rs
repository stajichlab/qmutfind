//! The locus drawing of an insertion: which genes it hits, what is lost, the copyable regions and the picture.

use qmutfind_core::{locus_view, parse_breakpoint_table, region_dots, region_text, Annotation, LocusInput, Orientation, Target};

/// A plus-strand gene G1 at 1,000-3,999 with two exons and a minus-strand gene G2 at 6,000-6,999, and a gene G3 that is far away.
const GFF: &str = "\
##gff-version 3
scf1\tt\tgene\t1000\t3999\t.\t+\t.\tID=G1;locus_tag=G1;Name=pfk
scf1\tt\tmRNA\t1000\t3999\t.\t+\t.\tID=G1-T;Parent=G1;product=phosphofructokinase
scf1\tt\texon\t1000\t1999\t.\t+\t.\tID=e1;Parent=G1-T
scf1\tt\texon\t2500\t3999\t.\t+\t.\tID=e2;Parent=G1-T
scf1\tt\tCDS\t1100\t1999\t.\t+\t0\tID=c1;Parent=G1-T
scf1\tt\tCDS\t2500\t3900\t.\t+\t0\tID=c2;Parent=G1-T
scf1\tt\tgene\t6000\t6999\t.\t-\t.\tID=G2
scf1\tt\tmRNA\t6000\t6999\t.\t-\t.\tID=G2-T;Parent=G2
scf1\tt\texon\t6000\t6999\t.\t-\t.\tID=e3;Parent=G2-T
scf1\tt\tCDS\t6050\t6950\t.\t-\t0\tID=c3;Parent=G2-T
scf1\tt\tgene\t90000\t90500\t.\t+\t.\tID=G3
scf1\tt\tmRNA\t90000\t90500\t.\t+\t.\tID=G3-T;Parent=G3
scf1\tt\texon\t90000\t90500\t.\t+\t.\tID=e4;Parent=G3-T
scf1\tt\tCDS\t90000\t90500\t.\t+\t0\tID=c4;Parent=G3-T
";

fn ann() -> Annotation {
    Annotation::from_gff(GFF.as_bytes()).unwrap()
}

fn input(left: u64, right: u64, inserted: u64) -> LocusInput {
    LocusInput {
        contig: "scf1".into(),
        left: Some(left),
        right: Some(right),
        inserted_bases: Some(inserted),
        label: "hph".into(),
        orientation: Some(Orientation::Forward),
        target: None,
        shared_left: 0,
        shared_right: 0,
    }
}

#[test]
fn removing_the_3_prime_end_of_a_plus_strand_gene_is_named_with_its_stop_codon() {
    // 3,001-4,200 is removed: the end of G1, with its stop codon at 3,900.
    let l = locus_view(&input(3000, 4201, 800), &ann(), None);
    assert_eq!(l.genes.len(), 1);
    let g = &l.genes[0];
    assert_eq!(g.id, "G1");
    assert_eq!(g.removed_bases, 999);
    assert_eq!(g.coding_bases, 900 + 1401);
    assert_eq!(g.coding_removed, 3900 - 3001 + 1);
    assert_eq!(g.effect, "3′ end removed; stop codon removed");
}

#[test]
fn the_same_loss_on_a_minus_strand_gene_is_its_5_prime_end_and_its_start_codon() {
    // The left end of a minus-strand gene is its 3' end: 5,500-6,500 is removed, with the stop codon at 6,050.
    let three = locus_view(&input(5499, 6501, 100), &ann(), None);
    assert_eq!(three.genes[0].effect, "3′ end removed; stop codon removed");
    // The right end of a minus gene is the 5' end, and the start codon sits at 6,950.
    let five = locus_view(&input(6499, 7501, 100), &ann(), None);
    assert_eq!(five.genes[0].effect, "5′ end removed; start codon removed");
}

#[test]
fn a_span_inside_a_gene_interrupts_it_and_a_whole_gene_is_named() {
    let inside = locus_view(&input(2099, 2301, 50), &ann(), None);
    assert_eq!(inside.genes[0].effect, "internal part removed, gene interrupted");
    let whole = locus_view(&input(900, 4100, 50), &ann(), None);
    assert_eq!(whole.genes[0].effect, "whole gene removed");
    assert_eq!(whole.genes[0].removed_bases, 3000);
}

#[test]
fn a_replacement_between_genes_hits_none_and_a_plain_insertion_in_a_gene_interrupts_it() {
    assert!(locus_view(&input(4999, 5001, 500), &ann(), None).genes.is_empty());
    // A plain insertion: the right flank is the next base.
    let l = locus_view(&input(2200, 2201, 500), &ann(), None);
    assert_eq!(l.genes.len(), 1);
    assert_eq!(l.genes[0].removed_bases, 0);
    assert_eq!(l.genes[0].effect, "interrupted by the insertion");
}

#[test]
fn regions_are_given_as_copyable_text_in_both_forms() {
    let l = locus_view(&input(3000, 4201, 800), &ann(), None);
    let removed = l.regions.iter().find(|r| r.name == "removed reference bases").unwrap();
    assert_eq!(removed.text(), "scf1:3001-4200");
    assert_eq!(removed.dots(), "scf1:3001..4200");
    let breakpoints = l.regions.iter().find(|r| r.name.starts_with("breakpoints")).unwrap();
    assert_eq!(breakpoints.text(), "scf1:3000-4201");
    assert!(l.regions.iter().any(|r| r.name == "drawn window"));
    assert_eq!(l.genes[0].region(), "scf1:1000-3999");
    assert!(l.genes_tsv().contains("\tscf1:1000-3999\t"));
    assert!(l.regions_tsv().contains("removed reference bases\tscf1:3001-4200\tscf1:3001..4200"));
    assert_eq!(region_text("a", 1, 2), "a:1-2");
    assert_eq!(region_dots("a", 1, 2), "a:1..2");
}

#[test]
fn the_picture_draws_the_replacement_with_the_same_scale_on_both_rows_and_links_only_hit_genes() {
    let mut i = input(3000, 4201, 800);
    i.target = Some(Target { contig: "scf1".into(), start: 3001, end: 4200 });
    let l = locus_view(&i, &ann(), Some("https://fungidb.org/gene/"));
    assert!(l.svg.contains("class=\"removed-span\" data-start=\"3001\" data-end=\"4200\""));
    assert!(l.svg.contains("class=\"insertion\" data-bases=\"800\""));
    assert!(l.svg.contains("class=\"target\""));
    assert!(l.svg.contains("class=\"replacement\""));
    assert!(l.svg.contains("href=\"https://fungidb.org/gene/G1\""));
    assert!(!l.svg.contains("gene/G2"), "a gene that is not hit gets no link");
    assert!(l.svg.contains("the allele is 400 bp shorter than the reference"));
    assert!(l.svg.contains("left flank ends at 3,000") && l.svg.contains("right flank starts at 4,201"));
    // One pixel per base on both rows: the insertion is 800/1200 of the width of the removed span.
    let width = |class: &str| -> f64 {
        let at = l.svg.find(&format!("class=\"{class}\"")).unwrap();
        let rest = &l.svg[at..];
        let w = rest.split("width=\"").nth(1).unwrap().split('"').next().unwrap();
        w.parse().unwrap()
    };
    let ratio = width("insertion") / width("removed-span");
    assert!((ratio - 800.0 / 1200.0).abs() < 0.01, "ratio {ratio}");
}

#[test]
fn a_gene_id_with_odd_characters_is_encoded_in_the_link_and_text_is_escaped() {
    let gff = "scf1\tt\tgene\t10\t100\t.\t+\t.\tID=a b<c;locus_tag=a b<c;Name=x\nscf1\tt\tmRNA\t10\t100\t.\t+\t.\tID=m;Parent=a b<c\nscf1\tt\texon\t10\t100\t.\t+\t.\tID=e;Parent=m\nscf1\tt\tCDS\t10\t100\t.\t+\t0\tID=c;Parent=m\n";
    let a = Annotation::from_gff(gff.as_bytes()).unwrap();
    let mut i = input(49, 51, 5);
    i.label = "<script>".into();
    let l = locus_view(&i, &a, Some("https://x/g/"));
    assert!(!l.svg.contains("<script>"));
    assert!(l.svg.contains("&lt;script&gt;"));
    assert!(l.svg.contains("href=\"https://x/g/a%20b%3Cc\""));
}

#[test]
fn without_genes_the_picture_still_shows_the_replacement() {
    let l = locus_view(&input(3000, 4201, 800), &Annotation::empty(), None);
    assert!(l.genes.is_empty());
    assert!(l.svg.contains("No gene models in this window."));
    assert!(l.svg.contains("class=\"insertion\""));
}

#[test]
fn the_json_carries_the_picture_the_genes_and_the_regions() {
    let l = locus_view(&input(3000, 4201, 800), &ann(), None);
    let json = l.to_json();
    assert!(json.starts_with("{\"svg\":\"<svg"));
    assert!(json.contains("\"genes\":[{\"id\":\"G1\",\"name\":\"pfk\",\"product\":\"phosphofructokinase\",\"region\":\"scf1:1000-3999\""));
    assert!(json.contains("\"region_dots\":\"scf1:3001..4200\""));
}

#[test]
fn a_breakpoint_table_is_read_by_column_name_and_gives_the_inputs() {
    let tsv = "sample\tcontig\tstrand\tref_scaffold\tref_left_end\tref_right_start\tref_gap_bp\tquery_gap_bp\tcall\toverlaps_locus\tmin_mapq\tinsert_fasta_id\n\
pfkA_control\tc9\t-\tscf1\t3000\t4201\t1200\t800\treplacement\tyes\t60\tx\n\
pfkA_control\tc9\t-\tscf1\t100\t170\t69\t66\treplacement\tno\t60\ty\n";
    let rows = parse_breakpoint_table(tsv).unwrap();
    assert_eq!(rows.len(), 2);
    assert!(rows[0].overlaps_locus && !rows[1].overlaps_locus);
    let i = rows[0].locus_input();
    assert_eq!((i.contig.as_str(), i.left, i.right, i.inserted_bases), ("scf1", Some(3000), Some(4201), Some(800)));
    assert_eq!(i.removed_span(), Some((3001, 4200)));
    // The gap that the table gives is the gap between the flanks.
    assert_eq!(rows[0].ref_gap_bp, rows[0].ref_right_start - rows[0].ref_left_end - 1);
    assert!(parse_breakpoint_table("sample\tcall\n").unwrap_err().contains("ref_scaffold"));
    assert!(parse_breakpoint_table("ref_scaffold\tref_left_end\tref_right_start\tref_gap_bp\tquery_gap_bp\nscf1\tx\t1\t1\t1\n").unwrap_err().contains("ref_left_end"));
}

#[test]
fn a_span_that_covers_exactly_the_coding_sequence_says_so_and_keeps_the_utr() {
    // G1 has CDS 1,100-1,999 and 2,500-3,900 and UTR at 1,000-1,099 and 3,901-3,999.
    let l = locus_view(&input(1099, 3901, 100), &ann(), None);
    let g = &l.genes[0];
    assert_eq!(g.coding_removed, g.coding_bases);
    assert_eq!(g.effect, "entire coding sequence removed (UTR sequence outside the span is kept)");
}

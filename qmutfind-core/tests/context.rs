//! The sequence around a variant: aligned rows of reference DNA, protein, variant DNA and protein, with exons in capitals.

mod common;

use common::{reverse_complement, Lcg};
use qmutfind_core::{sequence_context, Annotation, Genome, SequenceContext};

const AA: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

/// An independent translator for the expected values.
fn aa(codon: &[u8]) -> char {
    let i = |b: u8| match b { b'T' => 0, b'C' => 1, b'A' => 2, b'G' => 3, _ => panic!("base {b}") };
    AA[i(codon[0]) * 16 + i(codon[1]) * 4 + i(codon[2])] as char
}

fn protein(cds: &[u8]) -> Vec<char> {
    cds.chunks_exact(3).map(aa).collect()
}

/// A coding sequence of `len` bases (a multiple of 3): ATG, random codons without a stop, TAA.
fn cds_of(len: usize, seed: u64) -> Vec<u8> {
    let mut rng = Lcg(seed);
    let mut cds = b"ATG".to_vec();
    while cds.len() < len - 3 {
        let codon: Vec<u8> = (0..3).map(|_| rng.base()).collect();
        if !matches!(&codon[..], b"TAA" | b"TAG" | b"TGA") { cds.extend(codon); }
    }
    cds.extend(b"TAA");
    cds
}

struct Case { genome: Genome, ann: Annotation, seq: Vec<u8>, cds: Vec<u8> }

/// Contig `c1` of 300 bases of random background. A plus-strand gene: exon 1 is 41..=`e1_end`, the intron follows for 20 bases, exon 2 ends at 126
/// (its start makes the CDS 66 bases long). Without `with_gene` there is no gene.
fn plus_gene(e1_end: u64) -> Case {
    let mut rng = Lcg(7);
    let mut seq: Vec<u8> = (0..300).map(|_| rng.base()).collect();
    let cds = cds_of(66, 11);
    let n1 = (e1_end - 40) as usize;
    seq[40..40 + n1].copy_from_slice(&cds[..n1]);
    let e2_start = e1_end + 21; // intron e1_end+1 ..= e1_end+20
    seq[(e2_start - 1) as usize..126].copy_from_slice(&cds[n1..]);
    let gff = format!(
        "##gff-version 3\nc1\ts\tgene\t41\t126\t.\t+\t.\tID=g1;locus_tag=G1\nc1\ts\tmRNA\t41\t126\t.\t+\t.\tID=t1;Parent=g1\n\
         c1\ts\texon\t41\t{e1_end}\t.\t+\t.\tParent=t1\nc1\ts\texon\t{e2_start}\t126\t.\t+\t.\tParent=t1\n\
         c1\ts\tCDS\t41\t{e1_end}\t.\t+\t0\tParent=t1\nc1\ts\tCDS\t{e2_start}\t126\t.\t+\t.\tParent=t1\n"
    );
    Case { genome: Genome::new(vec![("c1".into(), seq.clone())]), ann: Annotation::from_gff(gff.as_bytes()).unwrap(), seq, cds }
}

/// The same kind of gene on the minus strand: exon 1 of the transcript is 290 down to 261 (30 bases), the intron 260..=241, exon 2 is 240 down to 205.
fn minus_gene() -> Case {
    let mut rng = Lcg(9);
    let mut seq: Vec<u8> = (0..300).map(|_| rng.base()).collect();
    let cds = cds_of(66, 13);
    let comp = |b: u8| reverse_complement(&[b])[0];
    for (i, &b) in cds[..30].iter().enumerate() { seq[(290 - i) - 1] = comp(b); }
    for (i, &b) in cds[30..].iter().enumerate() { seq[(240 - i) - 1] = comp(b); }
    let gff = "##gff-version 3\nc1\ts\tgene\t205\t290\t.\t-\t.\tID=g1;locus_tag=G1\nc1\ts\tmRNA\t205\t290\t.\t-\t.\tID=t1;Parent=g1\n\
               c1\ts\texon\t261\t290\t.\t-\t.\tParent=t1\nc1\ts\texon\t205\t240\t.\t-\t.\tParent=t1\n\
               c1\ts\tCDS\t261\t290\t.\t-\t0\tParent=t1\nc1\ts\tCDS\t205\t240\t.\t-\t.\tParent=t1\n";
    Case { genome: Genome::new(vec![("c1".into(), seq.clone())]), ann: Annotation::from_gff(gff.as_bytes()).unwrap(), seq, cds }
}

fn base(c: &Case, pos: u64) -> u8 { c.seq[(pos - 1) as usize] }
fn ctx(c: &Case, pos: u64, reference: &str, alt: &str, flank: u64) -> SequenceContext {
    sequence_context(&c.genome, &c.ann, "c1", pos, reference, alt, flank).unwrap()
}
fn row<'a>(x: &'a SequenceContext, label: &str) -> &'a str {
    &x.rows.iter().find(|(l, _)| l == label).unwrap_or_else(|| panic!("no row {label}: {:?}", x.rows.iter().map(|r| &r.0).collect::<Vec<_>>())).1
}
fn alt_base(b: u8) -> String { (if b == b'G' { 'T' } else { 'G' }).to_string() }
fn text(bytes: &[u8]) -> String { String::from_utf8(bytes.to_vec()).unwrap() }

#[test]
fn a_snv_in_an_exon_gives_aligned_rows_with_the_marker_under_the_changed_column() {
    let c = plus_gene(70);
    let r = base(&c, 50);
    let x = ctx(&c, 50, &(r as char).to_string(), &alt_base(r), 10); // window 40..=60: 21 columns
    for label in ["position", "marker", "reference", "ref protein", "variant", "var protein"] {
        assert_eq!(row(&x, label).len(), 21, "row {label} has one character per column");
    }
    assert_eq!(row(&x, "marker").find('^'), Some(10));
    assert_eq!(row(&x, "marker").matches('^').count(), 1);
    // The reference row: capitals in the exon (41..), a small letter outside it (40).
    let expect: String = (40..=60).map(|p| { let b = base(&c, p) as char; if p >= 41 { b } else { b.to_ascii_lowercase() } }).collect();
    assert_eq!(row(&x, "reference"), expect);
    assert_eq!(&row(&x, "variant")[..10], &row(&x, "reference")[..10]);
    assert_eq!(row(&x, "variant").as_bytes()[10] as char, alt_base(r).chars().next().unwrap());
    assert_eq!(&row(&x, "variant")[11..], &row(&x, "reference")[11..]);
}

#[test]
fn the_protein_letter_sits_under_the_middle_base_of_its_codon_and_the_variant_row_shows_the_new_residue() {
    let c = plus_gene(70);
    let r = base(&c, 50);
    let alt = alt_base(r);
    let x = ctx(&c, 50, &(r as char).to_string(), &alt, 10);
    let prot = protein(&c.cds);
    let ref_row: Vec<char> = row(&x, "ref protein").chars().collect();
    let var_row: Vec<char> = row(&x, "var protein").chars().collect();
    for p in 40u64..=60 {
        let col = (p - 40) as usize;
        let idx = p as i64 - 41; // coding index
        let middle = idx >= 0 && idx % 3 == 1;
        assert_eq!(ref_row[col], if middle { prot[(idx / 3) as usize] } else { ' ' }, "ref protein at {p}");
    }
    // The codon of position 50 is codon 3 (indices 9..12), middle base at position 51, column 11.
    let mut changed = c.cds[9..12].to_vec();
    changed[(50 - 41) % 3] = alt.as_bytes()[0];
    assert_eq!(var_row[11], aa(&changed));
    assert_ne!(var_row[11], ref_row[11], "the test needs a missense change; change the seed if this fails");
    for col in (0..21).filter(|&col| col != 11) { assert_eq!(var_row[col], ref_row[col], "only the changed codon differs (column {col})"); }
}

#[test]
fn a_gene_on_the_minus_strand_is_shown_on_the_gene_strand_so_the_protein_reads_left_to_right() {
    let c = minus_gene();
    let pos = 280;
    let r = base(&c, pos);
    let x = ctx(&c, pos, &(r as char).to_string(), &alt_base(r), 10); // genomic 270..=290, shown from 290 down to 270
    let shown: Vec<u8> = reverse_complement(&c.seq[269..290]);
    let expect: String = shown.iter().enumerate().map(|(i, &b)| { let p = 290 - i as u64; if p >= 261 { b as char } else { b.to_ascii_lowercase() as char } }).collect();
    assert_eq!(row(&x, "reference"), expect);
    assert_eq!(row(&x, "marker").find('^'), Some(10)); // 290 - 280
    let prot = protein(&c.cds);
    let ref_row: Vec<char> = row(&x, "ref protein").chars().collect();
    for p in (270..=290u64).rev() {
        let col = (290 - p) as usize;
        let idx = 290 - p; // coding index
        assert_eq!(ref_row[col], if idx % 3 == 1 { prot[(idx / 3) as usize] } else { ' ' }, "ref protein at {p}");
    }
    assert!(x.header.contains("minus") || x.header.contains("reverse"), "{}", x.header);
}

#[test]
fn an_in_frame_deletion_leaves_gap_columns_in_the_variant_row_and_marks_each_deleted_column() {
    let c = plus_gene(70);
    let reference = text(&c.seq[48..52]); // anchor 49, deleted 50..=52
    let alt = text(&c.seq[48..49]);
    let x = ctx(&c, 49, &reference, &alt, 8); // window 41..=60
    assert_eq!(row(&x, "marker").matches('^').count(), 3);
    assert_eq!(row(&x, "marker").find('^'), Some(9)); // pos 50 is column 9
    assert_eq!(&row(&x, "variant")[9..12], "---");
    assert_eq!(row(&x, "variant").len(), row(&x, "reference").len());
    // Deleting one codon: the variant protein is the reference protein with codon 3's letter (codon indices 9..12) missing.
    let alt_cds: Vec<u8> = [&c.cds[..9], &c.cds[12..]].concat();
    let prot = protein(&alt_cds);
    let var_row: Vec<char> = row(&x, "var protein").chars().collect();
    // After the deletion the coding index of position p is p - 41 - 3; its middle bases carry letters of the shortened protein.
    for p in 53u64..=60 {
        let idx = p as i64 - 41 - 3;
        let col = (p - 41) as usize;
        if idx % 3 == 1 { assert_eq!(var_row[col], prot[(idx / 3) as usize], "var protein at {p}"); }
    }
}

#[test]
fn a_one_base_insertion_adds_a_gap_column_to_the_reference_row_and_shifts_the_protein_after_it() {
    let c = plus_gene(70);
    let anchor = text(&c.seq[49..50]); // pos 50
    let alt = format!("{anchor}A");
    let x = ctx(&c, 50, &anchor, &alt, 10); // window 40..=60 plus one inserted column: 22
    assert_eq!(row(&x, "reference").len(), 22);
    assert_eq!(row(&x, "reference").as_bytes()[11], b'-', "the inserted base has a gap in the reference row");
    assert_eq!(row(&x, "variant").as_bytes()[11], b'A');
    assert_eq!(row(&x, "marker").find('^'), Some(11));
    // The altered CDS: the insertion goes after coding index 9 (position 50).
    let alt_cds: Vec<u8> = [&c.cds[..10], b"A", &c.cds[10..]].concat();
    let prot = protein(&alt_cds);
    let var_row: Vec<char> = row(&x, "var protein").chars().collect();
    // Variant columns after the insertion: column 12 is position 51 with alt coding index 11, and so on.
    for col in 12..22usize {
        let p = 51 + (col - 12) as u64;
        let idx = p as i64 - 41 + 1;
        if p <= 70 && idx % 3 == 1 { assert_eq!(var_row[col], prot[(idx / 3) as usize], "var protein at position {p}"); }
    }
}

#[test]
fn an_intron_variant_is_in_small_letters_and_the_protein_does_not_change() {
    let c = plus_gene(70); // intron 71..=90
    let r = base(&c, 80);
    let x = ctx(&c, 80, &(r as char).to_string(), &alt_base(r), 12); // 68..=92
    let ref_row = row(&x, "reference");
    assert_eq!(&ref_row[0..3], &text(&c.seq[67..70])); // 68..=70 are exon: capitals
    assert!(ref_row[3..23].chars().all(|ch| ch.is_ascii_lowercase()), "the intron 71..=90 is in small letters: {ref_row}");
    assert!(ref_row[23..].chars().all(|ch| ch.is_ascii_uppercase()), "91.. is exon: {ref_row}");
    assert_eq!(row(&x, "ref protein"), row(&x, "var protein"));
}

#[test]
fn a_codon_split_by_an_intron_is_translated_across_the_junction_and_shown_under_its_middle_base() {
    let c = plus_gene(71); // exon 1 is 41..=71 (31 bases), intron 72..=91, exon 2 is 92..=126: codon 10 is positions 71, 92, 93
    let r = base(&c, 71);
    let alt = alt_base(r);
    let x = ctx(&c, 71, &(r as char).to_string(), &alt, 40); // 31..=111, with the intron in the middle
    let ref_row: Vec<char> = row(&x, "ref protein").chars().collect();
    let var_row: Vec<char> = row(&x, "var protein").chars().collect();
    let col_of = |p: u64| (p - 31) as usize;
    assert_eq!(ref_row[col_of(92)], aa(&c.cds[30..33]), "the middle base of the split codon is position 92");
    let mut changed = c.cds[30..33].to_vec();
    changed[0] = alt.as_bytes()[0];
    assert_eq!(var_row[col_of(92)], aa(&changed));
    assert!(ref_row[col_of(72)..col_of(92)].iter().all(|&ch| ch == ' '), "no protein letters in the intron");
}

#[test]
fn outside_genes_the_sequence_is_in_small_letters_with_no_protein_rows() {
    let c = plus_gene(70);
    let r = base(&c, 200);
    let x = ctx(&c, 200, &(r as char).to_string(), &alt_base(r), 10);
    assert!(row(&x, "reference").chars().all(|ch| ch.is_ascii_lowercase()));
    assert!(x.rows.iter().all(|(l, _)| l != "ref protein" && l != "var protein"));
    assert!(x.notes.iter().any(|n| n.contains("not in a gene")), "{:?}", x.notes);
}

#[test]
fn the_window_is_cut_at_the_start_of_the_contig() {
    let c = plus_gene(70);
    let r = base(&c, 5);
    let x = ctx(&c, 5, &(r as char).to_string(), &alt_base(r), 40);
    assert_eq!(row(&x, "reference").len(), 5 + 40); // positions 1..=45
    assert_eq!(row(&x, "marker").find('^'), Some(4));
}

#[test]
fn the_text_and_the_html_hold_the_same_rows_and_the_html_marks_the_changed_column() {
    let c = plus_gene(70);
    let r = base(&c, 50);
    let x = ctx(&c, 50, &(r as char).to_string(), &alt_base(r), 10);
    let t = x.to_text();
    assert!(t.lines().any(|l| l.starts_with("reference") && l.ends_with(row(&x, "reference"))), "{t}");
    let h = x.to_html();
    assert!(h.contains("<span class=\"ctx-var\">"), "{h}");
    let stripped = {
        let mut out = String::new();
        let mut in_tag = false;
        for ch in h.chars() { match ch { '<' => in_tag = true, '>' => in_tag = false, c if !in_tag => out.push(c), _ => {} } }
        out.replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&")
    };
    assert_eq!(stripped.trim_end(), t.trim_end(), "the html is the text with spans");
}

#[test]
fn a_contig_that_is_not_in_the_reference_or_a_reference_base_that_does_not_match_is_an_error() {
    let c = plus_gene(70);
    assert!(sequence_context(&c.genome, &c.ann, "nope", 5, "A", "C", 10).is_err());
    let r = base(&c, 50);
    let wrong = if r == b'A' { "C" } else { "A" };
    let err = sequence_context(&c.genome, &c.ann, "c1", 50, wrong, "T", 10).unwrap_err();
    assert!(err.contains("reference"), "{err}");
}

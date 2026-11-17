//! Breakpoints of an assembly against the reference near a locus. The assemblies are simulated, so the right answer is known.

mod common;

use common::{reverse_complement, Lcg};
use qmutfind_core::{breakpoint_table_tsv, inserted_sequences_fasta, parse_breakpoint_table, BreakpointFinder, BreakpointOptions, Region};

fn dna(n: usize, seed: u64) -> Vec<u8> {
    let mut rng = Lcg(seed);
    (0..n).map(|_| rng.base()).collect()
}

fn fasta(name: &str, seq: &[u8]) -> Vec<u8> {
    format!(">{name}\n{}\n", String::from_utf8_lossy(seq)).into_bytes()
}

struct World {
    genome: Vec<u8>,
    insert: Vec<u8>,
    finder: BreakpointFinder,
}

fn world() -> World {
    let genome = dna(60_000, 1);
    let insert = dna(1_500, 2);
    let finder = BreakpointFinder::new(&fasta("chr", &genome), "asm5").unwrap();
    World { genome, insert, finder }
}

fn locus() -> Region {
    Region { contig: "chr".into(), start: 19_001, end: 23_000 }
}

/// A few bases at a junction can equal the reference by chance, which moves the breakpoint by that many bases.
fn near(value: u64, expected: u64) -> bool {
    value.abs_diff(expected) <= 8
}

#[test]
fn a_replacement_gives_both_flanks_the_sizes_and_the_inserted_sequence() {
    let w = world();
    // 20,001-22,500 (2,500 bases) are replaced by 1,500 other bases.
    let assembly = [&w.genome[..20_000], &w.insert[..], &w.genome[22_500..]].concat();

    let found = w.finder.assembly("mut", &fasta("ctg1", &assembly), &locus(), &BreakpointOptions::default()).unwrap();

    assert_eq!(found.contigs_near_locus, 1);
    assert_eq!(found.rows.len(), 1, "{:?}", found.rows);
    let r = &found.rows[0];
    assert_eq!((r.call.as_str(), r.ref_scaffold.as_str(), r.contig.as_str(), r.strand), ("replacement", "chr", "ctg1", '+'));
    assert!(r.overlaps_locus);
    assert!(near(r.ref_left_end, 20_000) && near(r.ref_right_start, 22_501), "{} {}", r.ref_left_end, r.ref_right_start);
    assert!(near(r.query_gap_bp, 1_500) && near(r.ref_gap_bp, 2_500));
    assert_eq!(r.ref_gap_bp, r.ref_right_start - r.ref_left_end - 1);
    let seq = r.inserted.as_ref().unwrap();
    assert_eq!(seq.len() as u64, r.query_gap_bp);
    let insert = String::from_utf8(w.insert.clone()).unwrap();
    assert!(insert.contains(&seq[3..seq.len() - 3]), "the inserted sequence is the simulated one");
}

#[test]
fn an_assembly_on_the_other_strand_gives_the_same_flanks_and_the_sequence_in_reference_orientation() {
    let w = world();
    let assembly = reverse_complement(&[&w.genome[..20_000], &w.insert[..], &w.genome[22_500..]].concat());

    let found = w.finder.assembly("rev", &fasta("ctg1", &assembly), &locus(), &BreakpointOptions::default()).unwrap();

    let r = found.rows.iter().find(|r| r.overlaps_locus).expect("a gap at the locus");
    assert_eq!(r.strand, '-');
    assert!(near(r.ref_left_end, 20_000) && near(r.ref_right_start, 22_501));
    let seq = r.inserted.as_ref().unwrap();
    assert!(String::from_utf8(w.insert.clone()).unwrap().contains(&seq[5..seq.len() - 5]));
}

#[test]
fn a_plain_insertion_and_a_plain_deletion_are_called_as_such() {
    let w = world();
    let insertion = [&w.genome[..20_000], &w.insert[..], &w.genome[20_000..]].concat();
    let found = w.finder.assembly("ins", &fasta("c", &insertion), &locus(), &BreakpointOptions::default()).unwrap();
    let r = found.rows.iter().find(|r| r.overlaps_locus).expect("a gap");
    assert_eq!(r.call, "insertion");
    assert!(r.ref_gap_bp <= 8 && near(r.query_gap_bp, 1_500), "{r:?}");

    let deletion = [&w.genome[..20_000], &w.genome[22_500..]].concat();
    let found = w.finder.assembly("del", &fasta("c", &deletion), &locus(), &BreakpointOptions::default()).unwrap();
    let r = found.rows.iter().find(|r| r.overlaps_locus).expect("a gap");
    assert_eq!(r.call, "deletion");
    assert!(r.query_gap_bp <= 8 && near(r.ref_gap_bp, 2_500), "{r:?}");
    assert!(r.inserted.is_none() || r.inserted.as_ref().unwrap().len() <= 8);
}

#[test]
fn an_unchanged_assembly_has_no_gap_and_an_assembly_elsewhere_is_not_near_the_locus() {
    let w = world();
    let same = w.finder.assembly("wt", &fasta("c", &w.genome), &locus(), &BreakpointOptions::default()).unwrap();
    assert!(same.rows.is_empty());
    assert_eq!(same.contigs_near_locus, 1);

    // Only the first 10 kb of the genome: it does not reach the locus.
    let short = w.finder.assembly("short", &fasta("c", &w.genome[..10_000]), &locus(), &BreakpointOptions { flank: 1_000, ..Default::default() }).unwrap();
    assert!(short.rows.is_empty());
    assert_eq!(short.contigs_near_locus, 0);
}

#[test]
fn the_table_and_the_fasta_are_written_and_the_table_reads_back() {
    let w = world();
    let assembly = [&w.genome[..20_000], &w.insert[..], &w.genome[22_500..]].concat();
    let found = w.finder.assembly("mut", &fasta("ctg1", &assembly), &locus(), &BreakpointOptions::default()).unwrap();

    let tsv = breakpoint_table_tsv(&found.rows);
    assert!(tsv.starts_with("sample\tcontig\tstrand\tref_scaffold\tref_left_end\tref_right_start\tref_gap_bp\tquery_gap_bp\tcall\toverlaps_locus\tmin_mapq\tinsert_fasta_id\n"));
    let back = parse_breakpoint_table(&tsv).unwrap();
    assert_eq!(back.len(), 1);
    assert_eq!((back[0].ref_left_end, back[0].ref_right_start, back[0].query_gap_bp), (found.rows[0].ref_left_end, found.rows[0].ref_right_start, found.rows[0].query_gap_bp));
    assert!(back[0].overlaps_locus && back[0].min_mapq == 60 && back[0].strand == '+');
    let fa = inserted_sequences_fasta(&found.rows);
    assert!(fa.starts_with(&format!(">{}\n", found.rows[0].insert_id)));
    assert!(found.rows[0].insert_id.starts_with("mut|ctg1:20") && found.rows[0].insert_id.contains("|ref_"));
}

#[test]
fn a_locus_that_is_not_in_the_reference_and_a_bad_preset_are_errors_that_say_why() {
    let w = world();
    let bad = Region { contig: "nope".into(), start: 1, end: 10 };
    let err = w.finder.assembly("x", &fasta("c", &w.genome), &bad, &BreakpointOptions::default()).unwrap_err();
    assert!(err.contains("nope"), "{err}");
    let err = BreakpointFinder::new(&fasta("chr", &w.genome), "asm7").err().unwrap();
    assert!(err.contains("asm5, asm10 or asm20"), "{err}");
    assert!(BreakpointFinder::new(b"", "asm5").err().unwrap().contains("no sequences"));
}

#[test]
fn a_plain_insertion_in_an_assembly_on_the_other_strand_is_cut_at_the_right_query_positions() {
    let w = world();
    let assembly = reverse_complement(&[&w.genome[..20_000], &w.insert[..], &w.genome[20_000..]].concat());
    let found = w.finder.assembly("rev", &fasta("c", &assembly), &locus(), &BreakpointOptions::default()).unwrap();
    let r = found.rows.iter().find(|r| r.overlaps_locus).expect("a gap");
    assert_eq!((r.call.as_str(), r.strand), ("insertion", '-'));
    assert!(near(r.query_gap_bp, 1_500) && near(r.ref_left_end, 20_000), "{r:?}");
    // The inserted bases are cut from the right place of the reverse strand contig and written as they read on the reference strand.
    let seq = r.inserted.as_ref().unwrap();
    assert!(String::from_utf8(w.insert.clone()).unwrap().contains(&seq[5..seq.len() - 5]));
}

#[test]
fn a_replacement_of_the_same_length_is_found_although_the_aligner_aligns_straight_through_it() {
    let w = world();
    // 2,000 bases are replaced by 2,000 other bases: no change in length, so the aligner has no reason to split the alignment.
    let other = dna(2_000, 9);
    let assembly = [&w.genome[..20_000], &other[..], &w.genome[22_000..]].concat();
    let found = w.finder.assembly("same", &fasta("c", &assembly), &locus(), &BreakpointOptions::default()).unwrap();
    let r = found.rows.iter().find(|r| r.overlaps_locus).unwrap_or_else(|| panic!("no gap: {:?}", found.rows));
    assert_eq!(r.call, "replacement");
    assert!(near(r.ref_left_end, 20_000) && near(r.ref_right_start, 22_001), "{r:?}");
    assert!(near(r.query_gap_bp, 2_000) && near(r.ref_gap_bp, 2_000), "{r:?}");
}

#[test]
fn a_few_snps_and_a_short_indel_are_not_gaps() {
    let w = world();
    let mut assembly = w.genome.clone();
    for i in (20_100..22_000).step_by(97) {
        assembly[i] = if assembly[i] == b'A' { b'C' } else { b'A' };
    }
    let assembly = [&assembly[..21_000], &dna(30, 4)[..], &assembly[21_000..]].concat();
    let found = w.finder.assembly("snps", &fasta("c", &assembly), &locus(), &BreakpointOptions::default()).unwrap();
    assert!(found.rows.is_empty(), "{:?}", found.rows);
}

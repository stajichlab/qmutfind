//! Behavior of the haploid pileup caller.

mod common;

use common::{other_base, reference};
use qmutfind_core::{Aligner, Caller, CallerConfig, ReadKind, RammapAligner};

/// Ten 150 bp reads that all cover 0-based position 50_000 of chrA, each carrying `alt` there.
fn reads_with_alt_at_50000(refs: &[(String, Vec<u8>)], alt: u8, count: usize) -> Vec<Vec<u8>> {
    (0..count)
        .map(|i| {
            let start = 49_900 + i * 10;
            let mut read = refs[0].1[start..start + 150].to_vec();
            read[50_000 - start] = alt;
            read
        })
        .collect()
}

fn config() -> CallerConfig {
    CallerConfig { min_depth: 5, min_alt_fraction: 0.5 }
}

#[test]
fn ten_reads_that_all_carry_a_snp_give_one_snp_call_with_depth_and_alt_count() {
    let refs = reference();
    let ref_base = refs[0].1[50_000];
    let alt = other_base(ref_base);
    let aligner = RammapAligner::new(refs.clone(), ReadKind::Short);
    let mut caller = Caller::new(refs.clone(), config());

    for (i, read) in reads_with_alt_at_50000(&refs, alt, 10).iter().enumerate() {
        for a in aligner.align(&format!("r{i}"), read) {
            caller.add(&a, read);
        }
    }
    let calls = caller.calls();

    assert_eq!(calls.len(), 1);
    let c = &calls[0];
    assert_eq!(c.contig, "chrA");
    assert_eq!(c.pos, 50_001); // 1-based
    assert_eq!(c.reference, (ref_base as char).to_string());
    assert_eq!(c.alt, (alt as char).to_string());
    assert_eq!(c.depth, 10);
    assert_eq!(c.alt_count, 10);
}

fn call_with(refs: &[(String, Vec<u8>)], reads: &[Vec<u8>], config: CallerConfig) -> Vec<qmutfind_core::Call> {
    let aligner = RammapAligner::new(refs.to_vec(), ReadKind::Short);
    let mut caller = Caller::new(refs.to_vec(), config);
    for (i, read) in reads.iter().enumerate() {
        for a in aligner.align(&format!("r{i}"), read) {
            caller.add(&a, read);
        }
    }
    caller.calls()
}

#[test]
fn an_alt_allele_in_only_two_of_ten_reads_is_not_called() {
    let refs = reference();
    let alt = other_base(refs[0].1[50_000]);
    let mut reads = reads_with_alt_at_50000(&refs, alt, 2);
    reads.extend(reads_with_alt_at_50000(&refs, refs[0].1[50_000], 8)); // 8 reads carry the reference base

    assert!(call_with(&refs, &reads, config()).is_empty());
}

#[test]
fn a_site_covered_by_fewer_reads_than_min_depth_is_not_called() {
    let refs = reference();
    let alt = other_base(refs[0].1[50_000]);
    let reads = reads_with_alt_at_50000(&refs, alt, 3); // all 3 carry alt, but min_depth is 5

    assert!(call_with(&refs, &reads, config()).is_empty());
}

/// Ten 150 bp reads that all carry a 2 bp deletion of 0-based positions 50_001 and 50_002.
fn reads_with_two_bp_deletion(refs: &[(String, Vec<u8>)], count: usize) -> Vec<Vec<u8>> {
    (0..count)
        .map(|i| {
            let start = 49_900 + i * 10;
            let mut read = refs[0].1[start..start + 152].to_vec();
            let at = 50_001 - start;
            read.drain(at..at + 2);
            read
        })
        .collect()
}

#[test]
fn ten_reads_with_a_two_bp_deletion_give_one_deletion_call_in_vcf_form() {
    let refs = reference();
    let reads = reads_with_two_bp_deletion(&refs, 10);

    let calls = call_with(&refs, &reads, config());

    assert_eq!(calls.len(), 1, "{calls:?}");
    let c = &calls[0];
    assert_eq!(c.contig, "chrA");
    assert_eq!(c.pos, 50_001); // 1-based position of the anchor base (0-based 50_000)
    let anchor_and_deleted = String::from_utf8(refs[0].1[50_000..50_003].to_vec()).unwrap();
    let anchor = String::from_utf8(refs[0].1[50_000..50_001].to_vec()).unwrap();
    assert_eq!(c.reference, anchor_and_deleted);
    assert_eq!(c.alt, anchor);
    assert_eq!(c.depth, 10);
    assert_eq!(c.alt_count, 10);
}

/// Ten 150 bp reads that all carry a 3 bp insertion between 0-based positions 50_000 and 50_001.
fn reads_with_three_bp_insertion(refs: &[(String, Vec<u8>)], inserted: &[u8], count: usize) -> Vec<Vec<u8>> {
    (0..count)
        .map(|i| {
            let start = 49_900 + i * 10;
            let mut read = refs[0].1[start..start + 147].to_vec();
            let at = 50_001 - start;
            for (k, &b) in inserted.iter().enumerate() {
                read.insert(at + k, b);
            }
            read
        })
        .collect()
}

#[test]
fn ten_reads_with_a_three_bp_insertion_give_one_insertion_call_in_vcf_form() {
    let refs = reference();
    let inserted = b"ACG";
    let reads = reads_with_three_bp_insertion(&refs, inserted, 10);

    let calls = call_with(&refs, &reads, config());

    assert_eq!(calls.len(), 1, "{calls:?}");
    let c = &calls[0];
    assert_eq!(c.contig, "chrA");
    assert_eq!(c.pos, 50_001);
    let anchor = String::from_utf8(refs[0].1[50_000..50_001].to_vec()).unwrap();
    assert_eq!(c.reference, anchor);
    assert_eq!(c.alt, format!("{anchor}ACG"));
    assert_eq!(c.depth, 10);
    assert_eq!(c.alt_count, 10);
}

#[test]
fn a_snp_at_a_site_where_the_reference_has_an_iupac_code_is_not_called() {
    let truth = reference();
    let mut ambiguous = truth.clone();
    ambiguous[0].1[50_000] = b'M'; // A or C; the real genome has such codes
    let reads = reads_with_alt_at_50000(&truth, other_base(truth[0].1[50_000]), 10);

    let calls = call_with(&ambiguous, &reads, config());

    assert!(calls.is_empty(), "{calls:?}");
}

#[test]
fn a_deletion_whose_reference_bases_include_an_iupac_code_is_not_called() {
    let truth = reference();
    let mut ambiguous = truth.clone();
    ambiguous[0].1[50_002] = b'R'; // inside the two deleted bases (50_001 and 50_002)
    let reads = reads_with_two_bp_deletion(&truth, 10);

    let calls = call_with(&ambiguous, &reads, config());

    assert!(calls.is_empty(), "{calls:?}");
}

/// Ten 3 kb reads that all carry a 100 bp insertion after 0-based position 71_500.
/// Each read has two substitution errors inside the inserted bases, at different places,
/// as real long reads do. Returns the reads and the true inserted sequence.
fn noisy_long_reads_with_insertion(refs: &[(String, Vec<u8>)]) -> (Vec<Vec<u8>>, Vec<u8>) {
    let mut rng = common::Lcg(7);
    let mut inserted: Vec<u8> = (0..100).map(|_| rng.base()).collect();
    // If the last inserted base equals the anchor base, the same insertion can also be written one
    // base to the left, and an aligner reports the leftmost form. Avoid that ambiguity here.
    let last = inserted.len() - 1;
    if inserted[last] == refs[0].1[71_500] {
        inserted[last] = other_base(inserted[last]);
    }
    let reads = (0..10)
        .map(|i| {
            let start = 70_000 + i * 20;
            let mut read = refs[0].1[start..start + 2_900].to_vec();
            let mut ins = inserted.clone();
            for pos in [(i * 7 + 3) % 100, (i * 13 + 50) % 100] {
                ins[pos] = other_base(ins[pos]);
            }
            let at = 71_501 - start;
            for (k, &b) in ins.iter().enumerate() {
                read.insert(at + k, b);
            }
            read
        })
        .collect();
    (reads, inserted)
}

#[test]
fn a_shared_insertion_is_called_with_the_consensus_sequence_despite_per_read_errors() {
    let refs = reference();
    let (reads, inserted) = noisy_long_reads_with_insertion(&refs);
    let aligner = RammapAligner::new(refs.clone(), ReadKind::Long);
    let mut caller = Caller::new(refs.clone(), config());

    for (i, read) in reads.iter().enumerate() {
        for a in aligner.align(&format!("l{i}"), read) {
            caller.add(&a, read);
        }
    }
    let calls = caller.calls();

    assert_eq!(calls.len(), 1, "{calls:?}");
    let c = &calls[0];
    assert_eq!(c.pos, 71_501);
    let anchor = (refs[0].1[71_500] as char).to_string();
    assert_eq!(c.reference, anchor);
    assert_eq!(c.alt, format!("{anchor}{}", String::from_utf8(inserted).unwrap()));
    assert_eq!(c.alt_count, 10);
}

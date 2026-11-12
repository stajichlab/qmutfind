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

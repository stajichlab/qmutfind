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
    CallerConfig { min_depth: 5, min_snp_fraction: 0.5, min_indel_fraction: 0.5, large_indel_min: 50 }
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
fn a_shared_large_insertion_is_called_as_a_symbolic_insertion_with_its_length() {
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
    assert_eq!(c.alt, "<INS>");
    assert_eq!(c.svlen, Some(inserted.len()));
    assert_eq!(c.alt_count, 10);
}

/// Ten 3 kb reads with an 8 bp insertion after 0-based position 71_500, one substitution error
/// per read inside the inserted bases, at different places.
fn noisy_long_reads_with_small_insertion(refs: &[(String, Vec<u8>)]) -> (Vec<Vec<u8>>, Vec<u8>) {
    let mut rng = common::Lcg(11);
    let mut inserted: Vec<u8> = (0..8).map(|_| rng.base()).collect();
    if inserted[7] == refs[0].1[71_500] {
        inserted[7] = other_base(inserted[7]);
    }
    let reads = (0..10)
        .map(|i| {
            let start = 70_000 + i * 20;
            let mut read = refs[0].1[start..start + 2_900].to_vec();
            let mut ins = inserted.clone();
            let bad = (i * 3 + 1) % 7; // never the last base
            ins[bad] = other_base(ins[bad]);
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
fn a_shared_small_insertion_is_called_with_the_consensus_sequence_despite_per_read_errors() {
    let refs = reference();
    let (reads, inserted) = noisy_long_reads_with_small_insertion(&refs);
    let aligner = RammapAligner::new(refs.clone(), ReadKind::Long);
    let mut caller = Caller::new(refs.clone(), config());

    for (i, read) in reads.iter().enumerate() {
        for a in aligner.align(&format!("s{i}"), read) {
            caller.add(&a, read);
        }
    }
    let calls = caller.calls();

    assert_eq!(calls.len(), 1, "{calls:?}");
    let anchor = (refs[0].1[71_500] as char).to_string();
    assert_eq!(calls[0].alt, format!("{anchor}{}", String::from_utf8(inserted).unwrap()));
    assert_eq!(calls[0].svlen, None);
}

/// Ten 3 kb reads that share one large insertion, but each read shows a slightly different
/// inserted length (94 to 106 bp), as real long reads do.
fn long_reads_with_jittered_insertion_lengths(refs: &[(String, Vec<u8>)]) -> Vec<Vec<u8>> {
    let mut rng = common::Lcg(5);
    let shared: Vec<u8> = (0..106).map(|_| rng.base()).collect();
    let lengths = [94usize, 96, 98, 100, 102, 104, 106, 95, 105, 100];
    lengths
        .iter()
        .enumerate()
        .map(|(i, &len)| {
            let start = 70_000 + i * 20;
            let mut read = refs[0].1[start..start + 2_900].to_vec();
            let at = 71_501 - start;
            for (k, &b) in shared[..len].iter().enumerate() {
                read.insert(at + k, b);
            }
            read
        })
        .collect()
}

#[test]
fn reads_with_slightly_different_insertion_lengths_support_one_symbolic_call_with_the_median_length() {
    let refs = reference();
    let reads = long_reads_with_jittered_insertion_lengths(&refs);
    let aligner = RammapAligner::new(refs.clone(), ReadKind::Long);
    let mut caller = Caller::new(refs.clone(), config());

    for (i, read) in reads.iter().enumerate() {
        for a in aligner.align(&format!("j{i}"), read) {
            caller.add(&a, read);
        }
    }
    let calls = caller.calls();

    assert_eq!(calls.len(), 1, "{calls:?}");
    let c = &calls[0];
    assert_eq!(c.alt, "<INS>");
    assert_eq!(c.svlen, Some(100));
    assert_eq!(c.alt_count, 10);
    // An aligner can place the same insertion a base or two apart, so allow a small window.
    assert!((71_497..=71_501).contains(&c.pos), "pos {}", c.pos);
}


fn config_with(snp: f64, indel: f64) -> CallerConfig {
    CallerConfig { min_depth: 5, min_snp_fraction: snp, min_indel_fraction: indel, large_indel_min: 50 }
}

#[test]
fn the_snp_threshold_applies_to_snps_and_the_indel_threshold_to_indels() {
    let refs = reference();
    let ref_base = refs[0].1[50_000];
    let alt = other_base(ref_base);
    // 7 of 10 reads carry the SNP, and 7 of 10 reads carry the 2 bp deletion: both at fraction 0.7.
    let mut snp_reads = reads_with_alt_at_50000(&refs, alt, 7);
    snp_reads.extend(reads_with_alt_at_50000(&refs, ref_base, 3));
    let mut del_reads = reads_with_two_bp_deletion(&refs, 7);
    del_reads.extend(reads_with_alt_at_50000(&refs, ref_base, 3));

    // A strict SNP threshold and a lenient indel threshold: the deletion is called, the SNP is not.
    assert!(call_with(&refs, &snp_reads, config_with(0.8, 0.5)).is_empty());
    assert_eq!(call_with(&refs, &del_reads, config_with(0.8, 0.5)).len(), 1);
    // The reverse: the SNP is called, the deletion is not.
    assert_eq!(call_with(&refs, &snp_reads, config_with(0.5, 0.8)).len(), 1);
    assert!(call_with(&refs, &del_reads, config_with(0.5, 0.8)).is_empty());
}

#[test]
fn the_defaults_depend_on_the_read_kind_and_are_stricter_for_snps_than_for_indels() {
    let short = CallerConfig::defaults_for(ReadKind::Short);
    let long = CallerConfig::defaults_for(ReadKind::Long);

    assert_eq!((short.min_depth, short.min_snp_fraction, short.min_indel_fraction, short.large_indel_min), (5, 0.8, 0.5, 50));
    assert_eq!((long.min_depth, long.min_snp_fraction, long.min_indel_fraction, long.large_indel_min), (5, 0.7, 0.5, 50));
}

#[test]
fn the_minimum_alt_reads_for_indels_applies_to_indels_but_not_to_snps() {
    let refs = reference();
    let ref_base = refs[0].1[50_000];
    let alt = other_base(ref_base);
    // 7 of 10 reads carry the SNP and 7 of 10 carry the deletion: seven supporting reads each.
    let mut snp_reads = reads_with_alt_at_50000(&refs, alt, 7);
    snp_reads.extend(reads_with_alt_at_50000(&refs, ref_base, 3));
    let mut del_reads = reads_with_two_bp_deletion(&refs, 7);
    del_reads.extend(reads_with_alt_at_50000(&refs, ref_base, 3));
    let run = |reads: &[Vec<u8>], min_alt_reads: u32| {
        let aligner = RammapAligner::new(refs.clone(), ReadKind::Short);
        let mut caller = Caller::new(refs.clone(), config_with(0.5, 0.5));
        caller.set_min_indel_alt_reads(min_alt_reads);
        for (i, read) in reads.iter().enumerate() {
            for a in aligner.align(&format!("r{i}"), read) {
                caller.add(&a, read);
            }
        }
        caller.calls().len()
    };

    assert_eq!(run(&del_reads, 7), 1); // seven reads is enough
    assert_eq!(run(&del_reads, 8), 0); // eight is not
    assert_eq!(run(&snp_reads, 8), 1); // the SNP is not an indel, so the setting does not apply
}

#[test]
fn the_default_minimum_alt_reads_for_indels_is_6_for_long_reads_and_0_for_short_reads() {
    assert_eq!(CallerConfig::default_min_indel_alt_reads(ReadKind::Long), 6);
    assert_eq!(CallerConfig::default_min_indel_alt_reads(ReadKind::Short), 0);
}

//! The confidence of a call: QUAL, a Phred-scaled binomial tail of the sample's sequencing-error rate, and
//! DPR, the depth at the call divided by the median depth of the sample (a repeat or duplication signal).

use qmutfind_core::{Call, CigarOp, Caller, CallerConfig, SamRecord};

const LEN: usize = 1000;

fn reference() -> Vec<(String, Vec<u8>)> {
    vec![("c".to_string(), (0..LEN).map(|i| b"ACGT"[(i * 7 + i / 3) % 4]).collect())]
}

fn other(base: u8) -> u8 {
    if base == b'A' { b'C' } else { b'A' }
}

fn config() -> CallerConfig {
    CallerConfig { min_depth: 1, min_snp_fraction: 0.5, min_indel_fraction: 0.5, large_indel_min: 50 }
}

/// A read that copies the reference from `start` for `len` bases, with the bases at `snps` (0-based) changed.
fn read(refs: &[(String, Vec<u8>)], start: usize, len: usize, snps: &[usize]) -> SamRecord {
    let mut seq = refs[0].1[start..start + len].to_vec();
    for &p in snps {
        seq[p - start] = other(seq[p - start]);
    }
    SamRecord { contig: "c".to_string(), ref_start: start, cigar_ops: vec![CigarOp::Match(len as u32)], seq, mapq: 60 }
}

fn calls_of(reads: Vec<SamRecord>) -> Vec<Call> {
    let mut caller = Caller::new(reference(), config());
    for r in &reads {
        caller.add_sam(r);
    }
    caller.calls()
}

/// P(X >= k) for X ~ Binomial(n, e), summed directly.
fn tail(n: u32, k: u32, e: f64) -> f64 {
    let choose = |n: u32, i: u32| (1..=i).fold(1.0, |acc, j| acc * (n - i + j) as f64 / j as f64);
    (k..=n).map(|i| choose(n, i) * e.powi(i as i32) * (1.0 - e).powi((n - i) as i32)).sum()
}

#[test]
fn qual_is_the_phred_scaled_binomial_tail_of_the_error_rate_estimated_from_the_data() {
    let refs = reference();
    // Five reads cover all 1000 bases. Four of them carry a change at 500. Of the 5,000 bases, 4 differ: e = 4/5000.
    let mut reads = vec![read(&refs, 0, LEN, &[])];
    reads.extend((0..4).map(|_| read(&refs, 0, LEN, &[500])));

    let calls = calls_of(reads);

    assert_eq!(calls.len(), 1);
    let expected = -10.0 * tail(5, 4, 4.0 / 5000.0).log10();
    let qual = calls[0].qual.expect("a pileup call has a QUAL");
    assert!((qual - expected).abs() < 0.05, "qual {qual}, expected {expected}");
}

#[test]
fn a_call_with_more_supporting_reads_has_a_higher_qual() {
    let refs = reference();
    let mut reads = Vec::new();
    for i in 0..5 {
        // Site 300: all five reads carry the change. Site 600: four of the five.
        let snps: Vec<usize> = if i < 4 { vec![300, 600] } else { vec![300] };
        reads.push(read(&refs, 0, LEN, &snps));
    }

    let calls = calls_of(reads);

    let q = |pos: u64| calls.iter().find(|c| c.pos == pos).unwrap().qual.unwrap();
    assert!(q(301) > q(601), "{} {}", q(301), q(601));
}

#[test]
fn the_depth_ratio_is_the_depth_at_the_call_over_the_median_depth() {
    let refs = reference();
    // Depth 5 on all 1000 bases (median 5). Five more reads cover 400 to 599 and carry a change at 500: depth 10 there.
    let mut reads: Vec<SamRecord> = (0..5).map(|_| read(&refs, 0, LEN, &[])).collect();
    reads.extend((0..5).map(|_| read(&refs, 400, 200, &[500])));

    let calls = calls_of(reads);

    assert_eq!(calls.len(), 1);
    assert_eq!(calls[0].depth, 10);
    let ratio = calls[0].depth_ratio.expect("a pileup call has a depth ratio");
    assert!((ratio - 2.0).abs() < 1e-9, "{ratio}");
}

#[test]
fn qual_is_finite_and_capped_when_the_tail_underflows() {
    let refs = reference();
    let reads: Vec<SamRecord> = (0..400).map(|_| read(&refs, 0, LEN, &[500])).collect();

    let calls = calls_of(reads);

    let qual = calls[0].qual.unwrap();
    assert!(qual.is_finite() && qual <= 999.0 && qual > 100.0, "{qual}");
}

#[test]
fn the_depth_ratio_uses_the_median_not_a_lower_quantile() {
    let refs = reference();
    // Depth 10 on 600 of the 1,000 bases and depth 5 on the other 400: the median is 10, so a call at depth 10 has a ratio of 1.
    let mut reads: Vec<SamRecord> = (0..5).map(|_| read(&refs, 0, LEN, &[])).collect();
    reads.extend((0..5).map(|_| read(&refs, 0, 600, &[500])));

    let calls = calls_of(reads);

    assert_eq!(calls.len(), 1);
    assert!((calls[0].depth_ratio.unwrap() - 1.0).abs() < 1e-9, "{:?}", calls[0].depth_ratio);
}

fn caller_with(reads: &[SamRecord]) -> Caller {
    let mut caller = Caller::new(reference(), config());
    for r in reads {
        caller.add_sam(r);
    }
    caller
}

#[test]
fn a_maximum_depth_ratio_drops_calls_in_regions_much_deeper_than_the_median() {
    let refs = reference();
    let mut reads: Vec<SamRecord> = (0..5).map(|_| read(&refs, 0, LEN, &[])).collect();
    reads.extend((0..5).map(|_| read(&refs, 400, 200, &[500]))); // depth 10 against a median of 5: ratio 2
    let mut caller = caller_with(&reads);

    assert_eq!(caller.calls().len(), 1, "no filter by default");
    caller.set_max_depth_ratio(1.9);
    assert_eq!(caller.calls().len(), 0);
    caller.set_max_depth_ratio(2.0);
    assert_eq!(caller.calls().len(), 1, "a ratio equal to the maximum is kept");
    assert_eq!(caller.max_depth_ratio(), 2.0);
}

#[test]
fn a_minimum_qual_drops_weakly_supported_calls_and_keeps_strong_ones() {
    let refs = reference();
    let mut reads = Vec::new();
    for i in 0..5 {
        let snps: Vec<usize> = if i < 3 { vec![300, 600] } else { vec![300] }; // 5 of 5 at 300, 3 of 5 at 600
        reads.push(read(&refs, 0, LEN, &snps));
    }
    let mut caller = caller_with(&reads);
    let calls = caller.calls();
    let (strong, weak) = (calls.iter().find(|c| c.pos == 301).unwrap().qual.unwrap(), calls.iter().find(|c| c.pos == 601).unwrap().qual.unwrap());
    assert!(strong > weak);

    caller.set_min_qual((strong + weak) / 2.0);

    let kept: Vec<u64> = caller.calls().iter().map(|c| c.pos).collect();
    assert_eq!(kept, vec![301]);
    assert_eq!(caller.min_qual(), (strong + weak) / 2.0);
}

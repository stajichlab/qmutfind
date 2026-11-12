//! SAM input: parse a record, and call variants from SAM records with the built-in caller.

mod common;

use common::{other_base, reference};
use qmutfind_core::{parse_sam_line, CigarOp, Caller, CallerConfig};

fn line(flag: u32, contig: &str, pos1: usize, cigar: &str, seq: &str) -> String {
    format!("r1\t{flag}\t{contig}\t{pos1}\t60\t{cigar}\t*\t0\t0\t{seq}\t{}", "I".repeat(seq.len()))
}

#[test]
fn a_record_gives_a_zero_based_start_and_the_cigar_operations() {
    let rec = parse_sam_line(&line(0, "chrA", 50_001, "4M2I3M1D2M", "ACGTTTACGAC")).unwrap().unwrap();

    assert_eq!(rec.contig, "chrA");
    assert_eq!(rec.ref_start, 50_000);
    assert_eq!(
        rec.cigar_ops,
        vec![CigarOp::Match(4), CigarOp::Insertion(2), CigarOp::Match(3), CigarOp::Deletion(1), CigarOp::Match(2)]
    );
    assert_eq!(rec.seq, b"ACGTTTACGAC".to_vec());
}

#[test]
fn clips_are_left_out_and_equals_and_x_count_as_matches() {
    // 2 soft-clipped bases, 3=, 1X, 2 more soft-clipped bases, 5 hard-clipped bases (not in SEQ).
    let rec = parse_sam_line(&line(16, "chrA", 11, "5H2S3=1X2S", "GGAAATCC")).unwrap().unwrap();

    assert_eq!(rec.cigar_ops, vec![CigarOp::Match(3), CigarOp::Match(1)]);
    assert_eq!(rec.seq, b"AAAT".to_vec()); // the aligned bases only, in reference orientation
}

#[test]
fn headers_unmapped_and_secondary_records_are_skipped() {
    assert_eq!(parse_sam_line("@SQ\tSN:chrA\tLN:1000"), Ok(None));
    assert_eq!(parse_sam_line(&line(4, "*", 0, "*", "ACGT")), Ok(None));       // 0x4 unmapped
    assert_eq!(parse_sam_line(&line(256, "chrA", 10, "4M", "ACGT")), Ok(None)); // 0x100 secondary
    assert!(parse_sam_line(&line(2048, "chrA", 10, "4M", "ACGT")).unwrap().is_some()); // supplementary is kept
}

#[test]
fn a_malformed_line_is_an_error_not_a_silent_skip() {
    assert!(parse_sam_line("r1\t0\tchrA\tnot-a-number\t60\t4M\t*\t0\t0\tACGT\tIIII").is_err());
    assert!(parse_sam_line("r1\t0\tchrA").is_err());
    assert!(parse_sam_line(&line(0, "chrA", 10, "4Q", "ACGT")).is_err());
}

#[test]
fn ten_sam_records_that_all_carry_a_snp_give_one_snp_call() {
    let refs = reference();
    let ref_base = refs[0].1[50_000];
    let alt = other_base(ref_base);
    let mut caller = Caller::new(
        refs.clone(),
        CallerConfig { min_depth: 5, min_alt_fraction: 0.5, large_indel_min: 50 },
    );

    for i in 0..10 {
        let start = 49_900 + i * 10;
        let mut seq = refs[0].1[start..start + 150].to_vec();
        seq[50_000 - start] = alt;
        let sam = line(0, "chrA", start + 1, "150M", &String::from_utf8(seq).unwrap());
        caller.add_sam(&parse_sam_line(&sam).unwrap().unwrap());
    }
    let calls = caller.calls();

    assert_eq!(calls.len(), 1);
    assert_eq!((calls[0].pos, calls[0].depth, calls[0].alt_count), (50_001, 10, 10));
    assert_eq!(calls[0].alt, (alt as char).to_string());
}

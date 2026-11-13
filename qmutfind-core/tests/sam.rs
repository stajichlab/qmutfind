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
        CallerConfig { min_depth: 5, min_snp_fraction: 0.5, min_indel_fraction: 0.5, large_indel_min: 50 },
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


fn line_with_mapq(mapq: u32, pos1: usize, seq: &str) -> String {
    format!("r1\t0\tchrA\t{pos1}\t{mapq}\t{}M\t*\t0\t0\t{seq}\t{}", seq.len(), "I".repeat(seq.len()))
}

#[test]
fn the_record_keeps_the_mapping_quality() {
    let rec = parse_sam_line(&line_with_mapq(37, 100, "ACGTAC")).unwrap().unwrap();

    assert_eq!(rec.mapq, 37);
}

#[test]
fn the_caller_ignores_sam_records_below_the_minimum_mapping_quality() {
    let refs = reference();
    let alt = other_base(refs[0].1[50_000]);
    let build = |min_mapq: u8| {
        let mut caller = Caller::new(
            refs.clone(),
            CallerConfig { min_depth: 5, min_snp_fraction: 0.5, min_indel_fraction: 0.5, large_indel_min: 50 },
        );
        caller.set_min_mapq(min_mapq);
        // Ten records that all carry the SNP: six with mapping quality 60, four with mapping quality 3.
        for i in 0..10 {
            let start = 49_900 + i * 10;
            let mut seq = refs[0].1[start..start + 150].to_vec();
            seq[50_000 - start] = alt;
            let mapq = if i < 6 { 60 } else { 3 };
            let sam = line_with_mapq(mapq, start + 1, &String::from_utf8(seq).unwrap());
            caller.add_sam(&parse_sam_line(&sam).unwrap().unwrap());
        }
        caller.calls()
    };

    assert_eq!(build(0)[0].depth, 10);
    assert_eq!(build(20)[0].depth, 6);
    assert!(build(61).is_empty());
}

/// A reference whose bases 50,000..50,006 (0-based) are a run of six A, between a C and a G.
fn reference_with_a_run() -> Vec<(String, Vec<u8>)> {
    let mut refs = reference();
    let seq = &mut refs[0].1;
    seq[49_999] = b'C';
    seq[50_000..50_006].copy_from_slice(b"AAAAAA");
    seq[50_006] = b'G';
    refs
}

fn config() -> CallerConfig {
    CallerConfig { min_depth: 5, min_snp_fraction: 0.5, min_indel_fraction: 0.5, large_indel_min: 50 }
}

#[test]
fn an_indel_in_a_homopolymer_run_is_called_in_its_leftmost_form_even_when_the_reads_place_it_at_the_right_end() {
    let refs = reference_with_a_run();
    let mut caller = Caller::new(refs.clone(), config());
    for i in 0..10 {
        // The reads delete the last A of the run (0-based 50,005), so their gap is at the right end of the run.
        let start = 49_950 - i;
        let mut seq = refs[0].1[start..50_005].to_vec();
        seq.extend_from_slice(&refs[0].1[50_006..50_060]);
        let sam = line(0, "chrA", start + 1, &format!("{}M1D54M", 50_005 - start), &String::from_utf8(seq).unwrap());
        caller.add_sam(&parse_sam_line(&sam).unwrap().unwrap());
    }

    let calls = caller.calls();

    assert_eq!(calls.len(), 1, "{calls:?}");
    // Leftmost form: anchor C at 1-based 50,000, then one A deleted.
    assert_eq!((calls[0].pos, calls[0].reference.as_str(), calls[0].alt.as_str()), (50_000, "CA", "C"));
}

#[test]
fn the_same_deletion_gives_the_same_call_whichever_base_of_the_run_the_reads_delete() {
    let refs = reference_with_a_run();
    let mut found = Vec::new();
    for deleted in [50_000usize, 50_003, 50_005] {
        let mut caller = Caller::new(refs.clone(), config());
        for i in 0..10 {
            let start = 49_950 - i;
            let mut seq = refs[0].1[start..deleted].to_vec();
            seq.extend_from_slice(&refs[0].1[deleted + 1..deleted + 55]);
            let sam = line(0, "chrA", start + 1, &format!("{}M1D54M", deleted - start), &String::from_utf8(seq).unwrap());
            caller.add_sam(&parse_sam_line(&sam).unwrap().unwrap());
        }
        let calls = caller.calls();
        assert_eq!(calls.len(), 1, "{calls:?}");
        found.push((calls[0].pos, calls[0].reference.clone(), calls[0].alt.clone()));
    }

    assert!(found.iter().all(|f| *f == found[0]), "{found:?}");
}

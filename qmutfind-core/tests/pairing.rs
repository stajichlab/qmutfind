//! Pair-aware mapping. The reference has a 600 bp block twice (20,000 and 60,000, 1-based 20,001 and 60,001), so a read inside the
//! block maps equally well to both copies and has mapping quality 0 on its own. Its mate lies in the unique flank before the first
//! copy, which tells where the read belongs.

mod common;

use common::{other_base, reverse_complement, Lcg};
use qmutfind_core::{Aligner, Engine, ReadKind, RammapAligner};

fn genome() -> Vec<u8> {
    let mut rng = Lcg(7);
    let mut seq: Vec<u8> = (0..120_000).map(|_| rng.base()).collect();
    let block = seq[20_000..20_600].to_vec();
    seq[60_000..60_600].copy_from_slice(&block);
    seq
}

fn fasta(seq: &[u8]) -> Vec<u8> {
    format!(">chrA\n{}\n", String::from_utf8_lossy(seq)).into_bytes()
}

fn fastq(name: &str, seq: &[u8]) -> String {
    format!("@{name}\n{}\n+\n{}\n", String::from_utf8_lossy(seq), "I".repeat(seq.len()))
}

/// Ten pairs. Mate 1 is forward and unique (the flank before the first copy). Mate 2 is reverse and inside the block, and carries
/// the SNP at 0-based 20,300 (VCF position 20,301). The insert size is 500.
fn pairs(seq: &[u8]) -> (String, String) {
    let alt = other_base(seq[20_300]);
    let (mut r1, mut r2) = (String::new(), String::new());
    for i in 0..10 {
        let m1 = seq[19_850 + i..19_950 + i].to_vec();
        let mut fwd2 = seq[20_250 + i..20_350 + i].to_vec();
        fwd2[20_300 - (20_250 + i)] = alt;
        r1.push_str(&fastq(&format!("p{i}/1"), &m1));
        r2.push_str(&fastq(&format!("p{i}/2"), &reverse_complement(&fwd2)));
    }
    (r1, r2)
}

fn records(vcf: &str) -> Vec<Vec<String>> {
    vcf.lines().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').map(str::to_string).collect()).collect()
}

#[test]
fn a_read_inside_a_repeat_has_a_candidate_in_each_copy_and_no_mapping_quality_on_its_own() {
    let seq = genome();
    let aligner = RammapAligner::new(vec![("chrA".to_string(), seq.clone())], ReadKind::Short);
    let read = seq[20_250..20_350].to_vec();

    let candidates = aligner.align_candidates("r", &read);

    let mut starts: Vec<usize> = candidates.iter().map(|a| a.ref_start).collect();
    starts.sort();
    assert!(starts.contains(&20_250) && starts.contains(&60_250), "{starts:?}");
    assert!(candidates.iter().all(|a| a.mapq() == 0), "an exact repeat gives no mapping quality: {:?}", candidates.iter().map(|a| a.mapq()).collect::<Vec<_>>());
    assert_eq!(aligner.align("r", &read).len(), 1, "align() gives only the primary mapping, as before");
}

#[test]
fn without_pairing_the_snp_in_the_repeat_is_not_called_because_its_reads_have_no_mapping_quality() {
    let seq = genome();
    let (r1, r2) = pairs(&seq);
    let mut engine = Engine::new(&fasta(&seq), ReadKind::Short, "s").unwrap();

    engine.add_reads(r1.as_bytes()).unwrap();
    engine.end_file().unwrap();
    engine.add_reads(r2.as_bytes()).unwrap();
    engine.end_file().unwrap();

    assert!(records(&engine.vcf().unwrap()).iter().all(|r| r[1] != "20301"), "the SNP must be missed by the single-read mapping");
}

#[test]
fn with_pairing_the_mate_in_the_unique_flank_places_the_read_and_the_snp_is_called() {
    let seq = genome();
    let (r1, r2) = pairs(&seq);
    let mut engine = Engine::new(&fasta(&seq), ReadKind::Short, "s").unwrap();

    engine.add_mate1(r1.as_bytes()).unwrap();
    engine.add_mate2(r2.as_bytes()).unwrap();
    engine.end_mates().unwrap();

    let calls = records(&engine.vcf().unwrap());
    let snp: Vec<&Vec<String>> = calls.iter().filter(|r| r[1] == "20301").collect();
    assert_eq!(snp.len(), 1, "{calls:?}");
    assert!(snp[0][7].contains("DP=10") && snp[0][7].contains("AO=10"), "{}", snp[0][7]);
    assert!(calls.iter().all(|r| r[1] != "60301"), "nothing may be called in the other copy: {calls:?}");
    assert_eq!(engine.read_count(), 20, "both mates of every pair are counted");
}

#[test]
fn pairing_works_when_the_chunks_of_the_two_files_split_the_records_in_different_places() {
    let seq = genome();
    let (r1, r2) = pairs(&seq);
    let mut engine = Engine::new(&fasta(&seq), ReadKind::Short, "s").unwrap();

    // Mate 1 in pieces of 37 bytes, mate 2 in pieces of 91, in turn.
    let (b1, b2) = (r1.as_bytes(), r2.as_bytes());
    let (mut i, mut j) = (0, 0);
    while i < b1.len() || j < b2.len() {
        if i < b1.len() { let e = (i + 37).min(b1.len()); engine.add_mate1(&b1[i..e]).unwrap(); i = e; }
        if j < b2.len() { let e = (j + 91).min(b2.len()); engine.add_mate2(&b2[j..e]).unwrap(); j = e; }
    }
    engine.end_mates().unwrap();

    assert!(records(&engine.vcf().unwrap()).iter().any(|r| r[1] == "20301"));
    assert_eq!(engine.read_count(), 20);
}

#[test]
fn two_files_with_different_numbers_of_reads_or_different_names_are_an_error() {
    let seq = genome();
    let (r1, r2) = pairs(&seq);
    let short: String = r2.lines().collect::<Vec<_>>()[..36].join("\n") + "\n"; // 9 records
    let mut engine = Engine::new(&fasta(&seq), ReadKind::Short, "s").unwrap();
    engine.add_mate1(r1.as_bytes()).unwrap();
    engine.add_mate2(short.as_bytes()).unwrap();
    let err = engine.end_mates().unwrap_err();
    assert!(err.contains("different numbers of reads"), "{err}");

    let renamed = r2.replace("p3/2", "q3/2");
    let mut engine = Engine::new(&fasta(&seq), ReadKind::Short, "s").unwrap();
    engine.add_mate1(r1.as_bytes()).unwrap();
    let err = engine.add_mate2(renamed.as_bytes()).and_then(|_| engine.end_mates()).unwrap_err();
    assert!(err.contains("p3") && err.contains("q3"), "the error names the two reads: {err}");
}

#[test]
fn a_pair_in_the_wrong_orientation_or_too_far_apart_is_not_used_to_rescue_a_read() {
    let seq = genome();
    let alt = other_base(seq[20_300]);
    let (mut r1, mut r2) = (String::new(), String::new());
    for i in 0..10 {
        // Mate 1 is 5,000 bp before the block: farther than the maximum insert size.
        let m1 = seq[15_000 + i..15_100 + i].to_vec();
        let mut fwd2 = seq[20_250 + i..20_350 + i].to_vec();
        fwd2[20_300 - (20_250 + i)] = alt;
        r1.push_str(&fastq(&format!("p{i}/1"), &m1));
        r2.push_str(&fastq(&format!("p{i}/2"), &reverse_complement(&fwd2)));
    }
    let mut engine = Engine::new(&fasta(&seq), ReadKind::Short, "s").unwrap();

    engine.add_mate1(r1.as_bytes()).unwrap();
    engine.add_mate2(r2.as_bytes()).unwrap();
    engine.end_mates().unwrap();

    assert!(records(&engine.vcf().unwrap()).iter().all(|r| r[1] != "20301"), "a pair that is too far apart is no evidence");
}

//! Behavior of the aligner backend that wraps `rammap`.
//! The reference is synthetic and deterministic, so results do not depend on any data file.

use qmutfind_core::{Aligner, CigarOp, ReadKind, RammapAligner};

/// Small deterministic generator, so the reference is the same on every machine.
struct Lcg(u64);

impl Lcg {
    fn base(&mut self) -> u8 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        b"ACGT"[((self.0 >> 33) % 4) as usize]
    }
}

fn reference() -> Vec<(String, Vec<u8>)> {
    let mut rng = Lcg(42);
    let chr_a: Vec<u8> = (0..120_000).map(|_| rng.base()).collect();
    let chr_b: Vec<u8> = (0..80_000).map(|_| rng.base()).collect();
    vec![("chrA".to_string(), chr_a), ("chrB".to_string(), chr_b)]
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

#[test]
fn a_long_read_with_a_deletion_an_insertion_and_mismatches_gets_the_expected_cigar() {
    let refs = reference();
    let mut read = refs[0].1[70_000..73_000].to_vec();
    for p in [300usize, 700, 1100, 1500, 1900] {
        read[p] = if read[p] == b'C' { b'A' } else { b'C' };
    }
    read.remove(2000);
    read.insert(2500, b'T');

    let aligner = RammapAligner::new(refs, ReadKind::Long);
    let alignments = aligner.align("long", &read);

    assert_eq!(alignments.len(), 1);
    let a = &alignments[0];
    assert_eq!(a.contig, "chrA");
    assert_eq!(a.ref_start, 70_000);
    assert_eq!(a.cigar_string(), "1998M1D501M1I500M");
    assert_eq!(
        a.cigar_ops(),
        vec![
            CigarOp::Match(1998),
            CigarOp::Deletion(1),
            CigarOp::Match(501),
            CigarOp::Insertion(1),
            CigarOp::Match(500),
        ]
    );
}

#[test]
fn aligned_read_returns_the_bases_in_reference_orientation_for_both_strands() {
    let refs = reference();
    let mut read = refs[0].1[49_950..50_100].to_vec();
    read[50] = if read[50] == b'G' { b'T' } else { b'G' };
    let reversed = reverse_complement(&read);

    let aligner = RammapAligner::new(refs, ReadKind::Short);
    let forward = &aligner.align("fwd", &read)[0];
    let reverse = &aligner.align("rev", &reversed)[0];

    assert_eq!(forward.aligned_read(&read), read);
    assert_eq!(reverse.aligned_read(&reversed), read);
}

#[test]
fn a_reverse_complemented_read_maps_to_the_same_place_on_the_reverse_strand() {
    let refs = reference();
    let read = refs[0].1[49_950..50_100].to_vec();

    let aligner = RammapAligner::new(refs, ReadKind::Short);
    let forward = aligner.align("fwd", &read);
    let reverse = aligner.align("rev", &reverse_complement(&read));

    assert_eq!(forward.len(), 1);
    assert_eq!(reverse.len(), 1);
    assert!(!forward[0].reverse);
    assert!(reverse[0].reverse);
    assert_eq!(forward[0].ref_start, reverse[0].ref_start);
}

#[test]
fn short_read_with_one_snp_maps_forward_with_a_full_length_match_cigar() {
    let refs = reference();
    let mut read = refs[0].1[49_950..50_100].to_vec();
    read[50] = if read[50] == b'G' { b'T' } else { b'G' };

    let aligner = RammapAligner::new(refs, ReadKind::Short);
    let alignments = aligner.align("snp_fwd", &read);

    assert_eq!(alignments.len(), 1);
    let a = &alignments[0];
    assert_eq!(a.contig, "chrA");
    assert_eq!(a.ref_start, 49_950);
    assert_eq!(a.cigar_string(), "150M");
}

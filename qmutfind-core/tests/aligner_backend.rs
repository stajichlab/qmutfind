//! Behavior of the aligner backend that wraps `rammap`.
//! The reference is synthetic and deterministic, so results do not depend on any data file.

use qmutfind_core::{Aligner, ReadKind, RammapAligner};

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

//! Helpers shared by the integration tests. The reference is synthetic and deterministic.
#![allow(dead_code)]

/// Small deterministic generator, so the reference is the same on every machine.
pub struct Lcg(pub u64);

impl Lcg {
    pub fn base(&mut self) -> u8 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        b"ACGT"[((self.0 >> 33) % 4) as usize]
    }
}

pub fn reference() -> Vec<(String, Vec<u8>)> {
    let mut rng = Lcg(42);
    let chr_a: Vec<u8> = (0..120_000).map(|_| rng.base()).collect();
    let chr_b: Vec<u8> = (0..80_000).map(|_| rng.base()).collect();
    vec![("chrA".to_string(), chr_a), ("chrB".to_string(), chr_b)]
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
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

/// A different base from `b`.
pub fn other_base(b: u8) -> u8 {
    if b == b'G' { b'T' } else { b'G' }
}

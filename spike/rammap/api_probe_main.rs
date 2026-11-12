use rammap::{Aligner, Preset, Strand};

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn base(&mut self) -> u8 { b"ACGT"[(self.next() % 4) as usize] }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', x => x }).collect()
}

fn show(tag: &str, name: &str, a: &Aligner, read: &[u8]) {
    let r = a.map_seq(name, read);
    if r.mappings.is_empty() { println!("{tag}\t{name}\tUNMAPPED"); return; }
    for m in &r.mappings {
        let s = if matches!(m.strand, Strand::Forward) { '+' } else { '-' };
        println!("{tag}\t{name}\t{}\t{}\t{}\t{}\t{}\t{}\tmapq={}\tprimary={}\tNM={}\tAS={}\tcigar={}\tmd={}",
            m.target_name, m.target_start, m.target_end, m.query_start, s, m.query_end, m.mapq,
            m.is_primary, m.edit_distance, m.score,
            m.cigar.as_deref().unwrap_or("-"), m.md.as_deref().unwrap_or("-"));
    }
}

fn main() {
    let mut rng = Lcg(42);
    let mut c1: Vec<u8> = (0..120_000).map(|_| rng.base()).collect();
    let c2: Vec<u8> = (0..80_000).map(|_| rng.base()).collect();
    c1[50_000] = b'A'; // keep a known base
    let seqs = vec![("chrA".to_string(), c1.clone()), ("chrB".to_string(), c2.clone())];

    // Illumina-like (sr): 150 bp reads with one SNP, and one with a 2 bp deletion and one with a 3 bp insertion.
    let sr = Aligner::from_seqs(seqs.clone(), Preset::Sr);
    let mut r1 = c1[49_950..50_100].to_vec();
    r1[50] = if r1[50] == b'G' { b'T' } else { b'G' };
    show("sr", "snp_fwd", &sr, &r1);
    show("sr", "snp_rev", &sr, &revcomp(&r1));
    let mut r2 = c2[10_000..10_160].to_vec();
    r2.drain(80..82);
    show("sr", "del2", &sr, &r2);
    let mut r3 = c2[30_000..30_150].to_vec();
    for (i, b) in b"ACG".iter().enumerate() { r3.insert(70 + i, *b); }
    show("sr", "ins3", &sr, &r3);

    // Nanopore-like (map-ont): 3 kb read with a 1 bp deletion, a 1 bp insertion and 5 mismatches.
    let ont = Aligner::from_seqs(seqs, Preset::MapOnt);
    let mut l = c1[70_000..73_000].to_vec();
    for p in [300usize, 700, 1100, 1500, 1900] { l[p] = if l[p] == b'C' { b'A' } else { b'C' }; }
    l.remove(2000);
    l.insert(2500, b'T');
    show("ont", "long_fwd", &ont, &l);
    show("ont", "long_rev", &ont, &revcomp(&l));
    // random read: should not map
    let junk: Vec<u8> = (0..2000).map(|_| rng.base()).collect();
    show("ont", "junk", &ont, &junk);
}

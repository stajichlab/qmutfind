//! Haploid pileup caller. Counts bases per reference position and calls sites where
//! one alternate allele has enough reads.

use std::collections::HashMap;

use crate::{Alignment, CigarOp};

/// Calling thresholds. Provisional defaults are `UNMEASURED` (see PLAN.md section 5).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CallerConfig {
    pub min_depth: u32,
    pub min_alt_fraction: f64,
}

/// One called variant. `pos` is 1-based.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Call {
    pub contig: String,
    pub pos: u64,
    pub reference: String,
    pub alt: String,
    pub depth: u32,
    pub alt_count: u32,
}

/// Reads that support one insertion, with base counts per inserted position.
struct InsertionEvidence {
    count: u32,
    bases: Vec<[u16; 4]>,
}

pub struct Caller {
    reference: Vec<(String, Vec<u8>)>,
    counts: Vec<Vec<[u16; 4]>>,
    /// (contig index, 0-based anchor position, deleted length) -> supporting reads
    deletions: HashMap<(usize, usize, usize), u32>,
    /// (contig index, 0-based anchor position, inserted length) -> evidence. Reads are grouped by
    /// length, not by exact sequence, because long reads carry errors inside the inserted bases.
    insertions: HashMap<(usize, usize, usize), InsertionEvidence>,
    config: CallerConfig,
}

fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

impl Caller {
    pub fn new(reference: Vec<(String, Vec<u8>)>, config: CallerConfig) -> Self {
        let counts = reference.iter().map(|(_, s)| vec![[0u16; 4]; s.len()]).collect();
        Caller { reference, counts, deletions: HashMap::new(), insertions: HashMap::new(), config }
    }

    /// Add one alignment. `read` is the sequence that was given to the aligner.
    pub fn add(&mut self, alignment: &Alignment, read: &[u8]) {
        let Some(ci) = self.reference.iter().position(|(n, _)| *n == alignment.contig) else {
            return;
        };
        let aligned = alignment.aligned_read(read);
        let (mut qi, mut ri) = (0usize, alignment.ref_start);
        for op in alignment.cigar_ops() {
            match *op {
                CigarOp::Match(n) => {
                    for i in 0..n as usize {
                        if let Some(b) = base_index(aligned[qi + i]) {
                            let slot = &mut self.counts[ci][ri + i][b];
                            *slot = slot.saturating_add(1);
                        }
                    }
                    qi += n as usize;
                    ri += n as usize;
                }
                CigarOp::Insertion(n) => {
                    if let Some(anchor) = ri.checked_sub(1) {
                        let ev = self
                            .insertions
                            .entry((ci, anchor, n as usize))
                            .or_insert_with(|| InsertionEvidence { count: 0, bases: vec![[0u16; 4]; n as usize] });
                        ev.count += 1;
                        for (i, &b) in aligned[qi..qi + n as usize].iter().enumerate() {
                            if let Some(bi) = base_index(b) {
                                ev.bases[i][bi] = ev.bases[i][bi].saturating_add(1);
                            }
                        }
                    }
                    qi += n as usize;
                }
                CigarOp::Deletion(n) => {
                    if let Some(anchor) = ri.checked_sub(1) {
                        *self.deletions.entry((ci, anchor, n as usize)).or_insert(0) += 1;
                    }
                    ri += n as usize;
                }
            }
        }
    }

    /// All calls, in reference order.
    pub fn calls(&self) -> Vec<Call> {
        let mut out: Vec<(usize, Call)> = Vec::new();
        for (ci, (name, seq)) in self.reference.iter().enumerate() {
            for (pos, c) in self.counts[ci].iter().enumerate() {
                let depth: u32 = c.iter().map(|&x| x as u32).sum();
                if depth < self.config.min_depth {
                    continue;
                }
                let Some(ref_idx) = base_index(seq[pos]) else {
                    continue; // no valid REF allele at an ambiguous reference base
                };
                let best = (0..4)
                    .filter(|&b| b != ref_idx)
                    .max_by_key(|&b| c[b])
                    .unwrap();
                let alt_count = c[best] as u32;
                if alt_count > 0 && alt_count as f64 / depth as f64 >= self.config.min_alt_fraction {
                    out.push((ci, Call {
                        contig: name.clone(),
                        pos: pos as u64 + 1,
                        reference: (seq[pos] as char).to_string(),
                        alt: (b"ACGT"[best] as char).to_string(),
                        depth,
                        alt_count,
                    }));
                }
            }
        }
        for (&(ci, anchor, len), &alt_count) in &self.deletions {
            let seq = &self.reference[ci].1;
            if !seq[anchor..anchor + 1 + len].iter().all(|&b| base_index(b).is_some()) {
                continue; // REF would hold a base that VCF does not allow
            }
            let depth: u32 = self.counts[ci][anchor].iter().map(|&x| x as u32).sum();
            if depth >= self.config.min_depth
                && alt_count as f64 / depth as f64 >= self.config.min_alt_fraction
            {
                out.push((ci, Call {
                    contig: self.reference[ci].0.clone(),
                    pos: anchor as u64 + 1,
                    reference: String::from_utf8_lossy(&seq[anchor..anchor + 1 + len]).into_owned(),
                    alt: (seq[anchor] as char).to_string(),
                    depth,
                    alt_count,
                }));
            }
        }
        for (&(ci, anchor, _len), ev) in &self.insertions {
            let alt_count = ev.count;
            let seq = &self.reference[ci].1;
            if base_index(seq[anchor]).is_none() {
                continue; // anchor base is ambiguous
            }
            let depth: u32 = self.counts[ci][anchor].iter().map(|&x| x as u32).sum();
            if depth >= self.config.min_depth
                && alt_count as f64 / depth as f64 >= self.config.min_alt_fraction
            {
                let mut alt = vec![seq[anchor]];
                alt.extend(ev.bases.iter().map(|c| {
                    let best = (0..4).max_by_key(|&b| c[b]).unwrap();
                    b"ACGT"[best]
                }));
                out.push((ci, Call {
                    contig: self.reference[ci].0.clone(),
                    pos: anchor as u64 + 1,
                    reference: (seq[anchor] as char).to_string(),
                    alt: String::from_utf8_lossy(&alt).into_owned(),
                    depth,
                    alt_count,
                }));
            }
        }
        out.sort_by_key(|(ci, c)| (*ci, c.pos));
        out.into_iter().map(|(_, c)| c).collect()
    }
}

/// Size class of a variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantClass {
    Snp,
    SmallIndel,
    LargeIndel,
}

/// Classify a variant given in VCF form. An indel has one anchor base plus the changed bases.
/// An indel is large when it changes `large_indel_min` bases or more.
pub fn classify(reference: &str, alt: &str, large_indel_min: usize) -> VariantClass {
    if reference.len() == alt.len() {
        return VariantClass::Snp;
    }
    let changed = reference.len().abs_diff(alt.len());
    if changed >= large_indel_min {
        VariantClass::LargeIndel
    } else {
        VariantClass::SmallIndel
    }
}

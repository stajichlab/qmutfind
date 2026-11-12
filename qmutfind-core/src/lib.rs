//! qmutfind core: pure Rust logic shared by the CLI and the WASM build.
//! This crate does no file I/O and starts no threads.

mod caller;
mod engine;
mod freebayes;
mod normalize;
mod parse;
mod sam;
mod vcf;
pub use caller::{classify, Call, Caller, CallerConfig, VariantClass};
pub use engine::Engine;
pub use freebayes::calls_from_freebayes_vcf;
pub use normalize::left_align_calls;
pub use parse::{read_fasta, FastqChunker, FastqReader};
pub use sam::{parse_sam_line, SamRecord};
pub use vcf::write_vcf;

use rammap::{Preset, Strand};

/// The kind of reads being mapped. It selects the aligner preset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadKind {
    Short,
    /// Nanopore-style long reads (the `map-ont` preset).
    Long,
}

/// One operation of an alignment. `Match` covers matches and mismatches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    Match(u32),
    Insertion(u32),
    Deletion(u32),
}

/// One alignment of a read to the reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub contig: String,
    /// 0-based start on the contig.
    pub ref_start: usize,
    /// True when the read aligned to the reverse strand.
    pub reverse: bool,
    cigar: String,
    cigar_ops: Vec<CigarOp>,
    query_start: usize,
    query_end: usize,
}

impl Alignment {
    pub fn cigar_string(&self) -> &str {
        &self.cigar
    }

    pub fn cigar_ops(&self) -> &[CigarOp] {
        &self.cigar_ops
    }

    /// The aligned part of `read`, in reference orientation. Clipped bases are left out.
    /// `read` is the sequence that was given to the aligner.
    pub fn aligned_read(&self, read: &[u8]) -> Vec<u8> {
        let part = &read[self.query_start..self.query_end];
        if self.reverse {
            part.iter().rev().map(|&b| complement(b)).collect()
        } else {
            part.to_vec()
        }
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// Maps reads to a reference. Backends implement this trait.
pub trait Aligner {
    fn align(&self, name: &str, seq: &[u8]) -> Vec<Alignment>;
}

/// Backend that uses the `rammap` library.
pub struct RammapAligner {
    inner: rammap::Aligner,
}

impl RammapAligner {
    pub fn new(reference: Vec<(String, Vec<u8>)>, kind: ReadKind) -> Self {
        let preset = match kind {
            ReadKind::Short => Preset::Sr,
            ReadKind::Long => Preset::MapOnt,
        };
        RammapAligner {
            inner: rammap::Aligner::from_seqs(reference, preset),
        }
    }
}

impl Aligner for RammapAligner {
    fn align(&self, name: &str, seq: &[u8]) -> Vec<Alignment> {
        self.inner
            .map_seq(name, seq)
            .mappings
            .iter()
            .filter(|m| m.is_primary)
            .map(|m| Alignment {
                contig: m.target_name.to_string(),
                ref_start: m.target_start,
                reverse: matches!(m.strand, Strand::Reverse),
                cigar: m.cigar.clone().unwrap_or_default(),
                cigar_ops: m
                    .cigar_ops
                    .iter()
                    .flatten()
                    .map(|c| match c.op {
                        0 => CigarOp::Match(c.len),
                        1 => CigarOp::Insertion(c.len),
                        2 => CigarOp::Deletion(c.len),
                        other => panic!(
                            "unsupported CIGAR op {other} from the aligner; DNA presets give only M, I and D"
                        ),
                    })
                    .collect(),
                query_start: m.query_start,
                query_end: m.query_end,
            })
            .collect()
    }
}

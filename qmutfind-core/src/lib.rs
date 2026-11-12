//! qmutfind core: pure Rust logic shared by the CLI and the WASM build.
//! This crate does no file I/O and starts no threads.

use rammap::Preset;

/// The kind of reads being mapped. It selects the aligner preset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadKind {
    Short,
}

/// One alignment of a read to the reference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    pub contig: String,
    /// 0-based start on the contig.
    pub ref_start: usize,
    cigar: String,
}

impl Alignment {
    pub fn cigar_string(&self) -> &str {
        &self.cigar
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
                cigar: m.cigar.clone().unwrap_or_default(),
            })
            .collect()
    }
}

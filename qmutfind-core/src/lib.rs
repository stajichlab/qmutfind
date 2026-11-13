//! qmutfind core: pure Rust logic shared by the CLI and the WASM build.
//! This crate does no file I/O and starts no threads.

mod annotate;
mod caller;
mod combine;
mod engine;
mod freebayes;
mod depth_plot;
mod genome_plot;
mod layout;
mod normalize;
mod pairing;
mod parse;
mod report;
mod sam;
mod stats;
mod vcf;
pub use combine::{combine, Cell, Combined, CombineInput, CombineOptions, Row, Status, StrainSummary};
pub use caller::{classify, Call, Caller, CallerConfig, DepthWindow, VariantClass};
pub use annotate::{annotate_text, annotate_vcf, classify_effect, Annotation, Effect, Genome, Transcript};
pub use engine::Engine;
pub use freebayes::calls_from_freebayes_vcf;
pub use depth_plot::depth_plot_svg;
pub use genome_plot::genome_plot_svg;
pub use normalize::left_align_calls;
pub use pairing::{choose_pair, DEFAULT_MAX_INSERT};
pub use parse::{read_fasta, FastqChunker, FastqReader};
pub use report::{
    parse_region, read_report_variants, render_locus_svg, render_report_html, report_page, Region, ReportVariant,
};
pub use sam::{parse_sam_line, SamRecord};
pub use stats::{compute_stats, ContigStats, Stats};
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
    mapq: u8,
    /// The alignment score of the aligner (higher is better).
    score: i32,
    /// False for a secondary alignment: another place that the read matches nearly as well as its primary one.
    primary: bool,
}

impl Alignment {
    pub fn cigar_string(&self) -> &str {
        &self.cigar
    }

    /// Mapping quality, 0 to 60 for rammap.
    pub fn mapq(&self) -> u8 {
        self.mapq
    }

    /// An alignment made by hand, for tests and for aligners other than `rammap`. `query_len` is the length of the read; the whole
    /// read is taken as aligned.
    pub fn new(contig: &str, ref_start: usize, reverse: bool, cigar_ops: Vec<CigarOp>, query_len: usize, mapq: u8, score: i32, primary: bool) -> Alignment {
        let cigar = cigar_ops
            .iter()
            .map(|op| match *op {
                CigarOp::Match(n) => format!("{n}M"),
                CigarOp::Insertion(n) => format!("{n}I"),
                CigarOp::Deletion(n) => format!("{n}D"),
            })
            .collect();
        Alignment { contig: contig.to_string(), ref_start, reverse, cigar, cigar_ops, query_start: 0, query_end: query_len, mapq, score, primary }
    }

    /// The alignment score of the aligner.
    pub fn score(&self) -> i32 {
        self.score
    }

    /// False for a secondary alignment.
    pub fn is_primary(&self) -> bool {
        self.primary
    }

    /// This alignment with another mapping quality (a pair can raise the quality of a read that is ambiguous on its own).
    pub fn with_mapq(mut self, mapq: u8) -> Alignment {
        self.mapq = mapq;
        self
    }

    /// The 0-based end on the contig (exclusive): the start plus the reference bases of the CIGAR.
    pub fn ref_end(&self) -> usize {
        self.ref_start
            + self
                .cigar_ops
                .iter()
                .map(|op| match *op {
                    CigarOp::Match(n) | CigarOp::Deletion(n) => n as usize,
                    CigarOp::Insertion(_) => 0,
                })
                .sum::<usize>()
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
    /// The primary alignments of a read (a read can have several when it is split).
    fn align(&self, name: &str, seq: &[u8]) -> Vec<Alignment>;

    /// All the places the read matches: the primary alignments and the secondary ones. Pairing chooses among them. A backend
    /// that has no secondary alignments gives the primary ones.
    fn align_candidates(&self, name: &str, seq: &[u8]) -> Vec<Alignment> {
        self.align(name, seq)
    }
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
        self.mappings(name, seq, false)
    }

    fn align_candidates(&self, name: &str, seq: &[u8]) -> Vec<Alignment> {
        self.mappings(name, seq, true)
    }
}

impl RammapAligner {
    fn mappings(&self, name: &str, seq: &[u8], with_secondary: bool) -> Vec<Alignment> {
        self.inner
            .map_seq(name, seq)
            .mappings
            .iter()
            .filter(|m| with_secondary || m.is_primary)
            .filter(|m| m.is_primary || m.cigar_ops.is_some()) // a secondary alignment without a CIGAR cannot be used
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
                mapq: m.mapq.clamp(0, 255) as u8,
                score: m.score,
                primary: m.is_primary,
            })
            .collect()
    }
}

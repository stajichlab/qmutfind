//! WebAssembly wrapper. A `Session` holds an `Engine`. JavaScript passes the reference FASTA
//! once, then FASTQ bytes in chunks, and reads the VCF at the end. Decompression (gzip, zstd)
//! is done in JavaScript before the bytes reach this code.

use qmutfind_core::{Engine, ReadKind};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct Session {
    engine: Engine,
}

fn read_kind(kind: &str) -> Result<ReadKind, JsError> {
    match kind {
        "short" => Ok(ReadKind::Short),
        "long" => Ok(ReadKind::Long),
        other => Err(JsError::new(&format!("kind must be short or long, not {other}"))),
    }
}

#[wasm_bindgen]
impl Session {
    /// `reference_fasta` is plain FASTA text as bytes. `kind` is "short" or "long".
    #[wasm_bindgen(constructor)]
    pub fn new(reference_fasta: &[u8], kind: &str, sample: &str) -> Result<Session, JsError> {
        let engine = Engine::new(reference_fasta, read_kind(kind)?, sample).map_err(|e| JsError::new(&e))?;
        Ok(Session { engine })
    }

    /// Replace the calling thresholds. Call before `add_reads`.
    pub fn set_thresholds(&mut self, min_depth: u32, min_snp_fraction: f64, min_indel_fraction: f64) {
        self.engine.set_thresholds(min_depth, min_snp_fraction, min_indel_fraction);
    }

    /// Ignore alignments with a mapping quality below `min_mapq`. Call before `add_reads`.
    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.engine.set_min_mapq(min_mapq);
    }

    /// Call an indel only when at least this many reads support it. Call before `add_reads`.
    pub fn set_min_indel_alt_reads(&mut self, min_alt_reads: u32) {
        self.engine.set_min_indel_alt_reads(min_alt_reads);
    }

    /// Drop calls with a QUAL below `min_qual`.
    pub fn set_min_qual(&mut self, min_qual: f64) {
        self.engine.set_min_qual(min_qual);
    }

    /// Drop calls deeper than `max_depth_ratio` times the median depth. 0 keeps all.
    pub fn set_max_depth_ratio(&mut self, max_depth_ratio: f64) {
        self.engine.set_max_depth_ratio(max_depth_ratio);
    }

    /// Add a chunk of FASTQ bytes. A chunk can split a record anywhere.
    pub fn add_reads(&mut self, chunk: &[u8]) -> Result<(), JsError> {
        self.engine.add_reads(chunk).map_err(|e| JsError::new(&e))
    }

    /// Mark the end of one reads file. A truncated last record is an error.
    pub fn end_file(&mut self) -> Result<(), JsError> {
        self.engine.end_file().map_err(|e| JsError::new(&e))
    }

    /// Reads mapped so far.
    pub fn read_count(&self) -> f64 {
        self.engine.read_count() as f64
    }

    /// The mean depth in windows of `window` bases as BED text with four columns (contig, start, end, depth).
    pub fn depth_bed(&self, window: u32) -> String {
        self.engine.depth_bed(window as usize)
    }

    /// The weak signals so far (calls at relaxed thresholds) as VCF text, for `combine_tables`.
    pub fn evidence_vcf(&self) -> Result<String, JsError> {
        self.engine.evidence_vcf().map_err(|e| JsError::new(&e))
    }

    /// The calls so far as VCF text.
    pub fn vcf(&self) -> Result<String, JsError> {
        self.engine.vcf().map_err(|e| JsError::new(&e))
    }
}

/// Add the effect of each variant to VCF text. `reference_fasta` and `gff` are plain bytes
/// (decompress in JavaScript first). Pass `undefined` for `gff` to mark every record unannotated.
#[wasm_bindgen]
pub fn annotate(vcf: &str, reference_fasta: &[u8], gff: Option<Vec<u8>>) -> Result<String, JsError> {
    qmutfind_core::annotate_text(vcf, reference_fasta, gff.as_deref()).map_err(|e| JsError::new(&e))
}

/// Confirmation of an introduced allele (see `qmutfind confirm`): the reads are mapped to the genome plus the construct, and the report says where
/// the construct joins the genome and whether it is intact. JavaScript passes the genome and the construct FASTA once, then FASTQ bytes in chunks.
#[wasm_bindgen]
pub struct Confirmation {
    inner: qmutfind_core::Confirmer,
}

#[wasm_bindgen]
impl Confirmation {
    /// `target_contig` is the contig of the region that the construct should replace (1-based `target_start` to `target_end`), or `undefined` for no target.
    #[wasm_bindgen(constructor)]
    pub fn new(
        genome_fasta: &[u8], construct_fasta: &[u8], kind: &str, target_contig: Option<String>, target_start: f64, target_end: f64,
        min_junction_reads: u32, min_depth: f64, tolerance: f64,
    ) -> Result<Confirmation, JsError> {
        let target = match target_contig {
            Some(contig) => {
                if !(target_start >= 1.0 && target_end >= target_start && target_end < 4.0e9) {
                    return Err(JsError::new("the target must be a region with 1 <= start <= end"));
                }
                Some(qmutfind_core::Target { contig, start: target_start as u64, end: target_end as u64 })
            }
            None => None,
        };
        let defaults = qmutfind_core::ConfirmOptions::default();
        let options = qmutfind_core::ConfirmOptions {
            target, min_junction_reads, min_depth, target_tolerance: tolerance.max(0.0) as u64, ..defaults
        };
        let inner = qmutfind_core::Confirmer::from_fasta(genome_fasta, construct_fasta, read_kind(kind)?, options).map_err(|e| JsError::new(&e))?;
        Ok(Confirmation { inner })
    }

    /// Add a chunk of FASTQ bytes. The chunk can split a record anywhere.
    pub fn add_reads(&mut self, chunk: &[u8]) -> Result<(), JsError> {
        self.inner.add_reads(chunk).map_err(|e| JsError::new(&e))
    }

    /// Mark the end of one reads file.
    pub fn end_file(&mut self) -> Result<(), JsError> {
        self.inner.end_file().map_err(|e| JsError::new(&e))
    }

    pub fn read_count(&self) -> f64 {
        self.inner.read_count() as f64
    }

    /// The report as text (the same as `PREFIX.confirm.txt` of the CLI).
    pub fn report_text(&self) -> String {
        self.inner.report().to_text()
    }

    /// The report as JSON.
    pub fn report_json(&self) -> String {
        self.inner.report().to_json()
    }

    /// All junctions as a table.
    pub fn junctions_tsv(&self) -> String {
        self.inner.report().junctions_tsv()
    }

    /// The depth at each base of the construct.
    pub fn depth_tsv(&self) -> String {
        self.inner.report().depth_tsv()
    }
}

/// The region of a gene in a GFF3 file (plain bytes), found by its ID, Name, gene or locus_tag, as JSON: `{"contig":..,"start":..,"end":..}`.
#[wasm_bindgen]
pub fn find_gene(gff: &[u8], name: &str) -> Result<String, JsError> {
    let t = qmutfind_core::find_gene(gff, name).map_err(|e| JsError::new(&e))?;
    Ok(format!("{{\"contig\":\"{}\",\"start\":{},\"end\":{}}}", t.contig.replace('\\', "\\\\").replace('"', "\\\""), t.start, t.end))
}

/// The sequence around one variant as aligned text rows (or as HTML with the changed columns in spans). `reference_fasta` and `gff` are plain
/// bytes; `reference` and `alt` are in VCF form. Same text as `qmutfind context`.
#[wasm_bindgen]
pub fn context(reference_fasta: &[u8], gff: Option<Vec<u8>>, contig: &str, pos: f64, reference: &str, alt: &str, flank: u32, html: bool) -> Result<String, JsError> {
    if !(pos >= 1.0 && pos.fract() == 0.0 && pos < 4.0e9) {
        return Err(JsError::new("the position must be a whole number of 1 or more"));
    }
    qmutfind_core::context_text(reference_fasta, gff.as_deref(), contig, pos as u64, reference, alt, flank as u64, html).map_err(|e| JsError::new(&e))
}

/// Summary statistics of VCF text as JSON. `window` is the distance in bases for the clustering count.
#[wasm_bindgen]
pub fn stats(vcf: &str, window: u32) -> Result<String, JsError> {
    qmutfind_core::compute_stats(vcf.as_bytes(), window as u64)
        .map(|s| s.to_json())
        .map_err(|e| JsError::new(&e))
}

/// The version of qmutfind (this package), for the records of a run.
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// `None` for an empty text.
fn non_empty(text: &str) -> Option<&str> {
    if text.is_empty() { None } else { Some(text) }
}

/// The combined table of several strains, its summary, the warnings (one per line, empty when there are none) and the merged VCF
/// (one record per distinct variant), as four texts [table, summary, warnings, vcf]. `names`, `vcfs` and `beds`
/// `evidence` have one entry per strain; an empty coverage or evidence text means the strain has none. `parent` is optional.
#[wasm_bindgen]
pub fn combine_tables(names: Vec<String>, vcfs: Vec<String>, beds: Vec<String>, evidence: Vec<String>, parent: Option<String>, min_depth: f64) -> Result<Vec<String>, JsError> {
    if names.len() != vcfs.len() || names.len() != beds.len() || names.len() != evidence.len() {
        return Err(JsError::new("give the same number of names, VCF texts, coverage texts and evidence texts"));
    }
    let inputs: Vec<qmutfind_core::CombineInput> = (0..names.len())
        .map(|i| qmutfind_core::CombineInput { name: &names[i], vcf: &vcfs[i], depth_bed: non_empty(&beds[i]), evidence_vcf: non_empty(&evidence[i]) })
        .collect();
    let combined = qmutfind_core::combine(&inputs, &qmutfind_core::CombineOptions { parent: parent.as_deref(), min_depth }).map_err(|e| JsError::new(&e))?;
    Ok(vec![combined.to_tsv(), combined.summary_tsv(), combined.warnings.join("\n"), combined.to_vcf()])
}

/// The locus-view page (HTML with an inline SVG) for one region, such as "c1:21-110".
/// Pass `undefined` for `gff` to draw the variants without gene models.
#[wasm_bindgen]
pub fn report(vcf: &str, gff: Option<Vec<u8>>, region: &str, title: &str) -> Result<String, JsError> {
    qmutfind_core::report_page(vcf, gff.as_deref(), region, title).map_err(|e| JsError::new(&e))
}

/// The whole-genome plot (SVG text) of VCF text. Contig lengths come from the VCF's `##contig` lines, or from
/// `reference_fasta` (plain FASTA bytes) when it is given.
#[wasm_bindgen]
pub fn genome_plot(vcf: &str, reference_fasta: Option<Vec<u8>>, order: Option<String>) -> Result<String, JsError> {
    let contigs = match reference_fasta {
        Some(bytes) => {
            let reference = qmutfind_core::read_fasta(bytes.as_slice()).map_err(|e| JsError::new(&format!("reference FASTA: {e}")))?;
            Some(reference.into_iter().map(|(name, seq)| (name, seq.len())).collect::<Vec<_>>())
        }
        None => None,
    };
    qmutfind_core::genome_plot_svg(vcf, contigs.as_deref(), order.as_deref().unwrap_or("input")).map_err(|e| JsError::new(&e))
}

/// The genome-wide depth plot (SVG text) from BED text (four columns; the output of `depth_bed` or of `mosdepth`).
/// `vcf` adds a tick for each call. `reference_fasta` (plain FASTA bytes) overrides the contig lengths.
/// `order` is "input", "length" or "name".
#[wasm_bindgen]
pub fn depth_plot(bed: &str, vcf: Option<String>, reference_fasta: Option<Vec<u8>>, order: Option<String>) -> Result<String, JsError> {
    let contigs = match reference_fasta {
        Some(bytes) => {
            let reference = qmutfind_core::read_fasta(bytes.as_slice()).map_err(|e| JsError::new(&format!("reference FASTA: {e}")))?;
            Some(reference.into_iter().map(|(name, seq)| (name, seq.len())).collect::<Vec<_>>())
        }
        None => None,
    };
    qmutfind_core::depth_plot_svg(bed, contigs.as_deref(), order.as_deref().unwrap_or("input"), vcf.as_deref()).map_err(|e| JsError::new(&e))
}

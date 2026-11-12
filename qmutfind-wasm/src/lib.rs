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

/// Summary statistics of VCF text as JSON. `window` is the distance in bases for the clustering count.
#[wasm_bindgen]
pub fn stats(vcf: &str, window: u32) -> Result<String, JsError> {
    qmutfind_core::compute_stats(vcf.as_bytes(), window as u64)
        .map(|s| s.to_json())
        .map_err(|e| JsError::new(&e))
}

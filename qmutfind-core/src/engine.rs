//! The whole pipeline over bytes: a reference FASTA and FASTQ chunks go in, VCF text comes out.
//! No file access and no threads. The browser build and the CLI both drive this type.

use crate::{
    read_fasta, write_vcf, Aligner, Caller, CallerConfig, FastqChunker, ReadKind, RammapAligner,
};

pub struct Engine {
    aligner: RammapAligner,
    caller: Caller,
    chunker: FastqChunker,
    contigs: Vec<(String, usize)>,
    sample: String,
    config: CallerConfig,
    reads: u64,
}

impl Engine {
    /// Build the index from `reference_fasta` (plain FASTA bytes). `kind` picks the aligner
    /// preset and the default thresholds.
    pub fn new(reference_fasta: &[u8], kind: ReadKind, sample: &str) -> Result<Engine, String> {
        let reference = read_fasta(reference_fasta).map_err(|e| format!("reference FASTA: {e}"))?;
        if reference.is_empty() {
            return Err("the reference FASTA has no sequences".to_string());
        }
        let contigs = reference.iter().map(|(n, s)| (n.clone(), s.len())).collect();
        let config = CallerConfig::defaults_for(kind);
        let mut caller = Caller::new(reference.clone(), config);
        caller.set_min_mapq(CallerConfig::default_min_mapq(kind));
        Ok(Engine {
            aligner: RammapAligner::new(reference, kind),
            caller,
            chunker: FastqChunker::new(),
            contigs,
            sample: sample.to_string(),
            config,
            reads: 0,
        })
    }

    /// Replace the calling thresholds. Call before adding reads.
    pub fn set_thresholds(&mut self, min_depth: u32, min_snp_fraction: f64, min_indel_fraction: f64) {
        self.config = CallerConfig { min_depth, min_snp_fraction, min_indel_fraction, ..self.config };
        self.caller.set_config(self.config);
    }

    /// Ignore alignments with a mapping quality below `min_mapq`. Call before adding reads.
    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.caller.set_min_mapq(min_mapq);
    }

    /// The minimum mapping quality in use.
    pub fn min_mapq(&self) -> u8 {
        self.caller.min_mapq()
    }

    /// Add a chunk of FASTQ bytes. The chunk can split a record anywhere.
    pub fn add_reads(&mut self, chunk: &[u8]) -> Result<(), String> {
        let records = self.chunker.push(chunk).map_err(|e| format!("reads: {e}"))?;
        self.map_records(records);
        Ok(())
    }

    /// Mark the end of one reads file. A truncated last record is an error.
    pub fn end_file(&mut self) -> Result<(), String> {
        let records = self.chunker.finish().map_err(|e| format!("reads: {e}"))?;
        self.map_records(records);
        Ok(())
    }

    fn map_records(&mut self, records: Vec<(String, Vec<u8>)>) {
        for (name, seq) in records {
            for alignment in self.aligner.align(&name, &seq) {
                self.caller.add(&alignment, &seq);
            }
            self.reads += 1;
        }
    }

    /// Reads mapped so far.
    pub fn read_count(&self) -> u64 {
        self.reads
    }

    /// The calls so far as VCF text.
    pub fn vcf(&self) -> Result<String, String> {
        let mut out = Vec::new();
        write_vcf(&mut out, &self.contigs, &self.sample, &self.caller.calls(), self.config.large_indel_min)
            .map_err(|e| e.to_string())?;
        String::from_utf8(out).map_err(|e| e.to_string())
    }
}

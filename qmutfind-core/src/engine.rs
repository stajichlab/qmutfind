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
        caller.set_min_indel_alt_reads(CallerConfig::default_min_indel_alt_reads(kind));
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

    /// Call an indel only when at least this many reads support it. Call before adding reads.
    pub fn set_min_indel_alt_reads(&mut self, min_alt_reads: u32) {
        self.caller.set_min_indel_alt_reads(min_alt_reads);
    }

    /// The minimum number of reads for an indel in use.
    pub fn min_indel_alt_reads(&self) -> u32 {
        self.caller.min_indel_alt_reads()
    }

    /// Drop calls with a QUAL below `min_qual`. Call before reading the VCF.
    pub fn set_min_qual(&mut self, min_qual: f64) {
        self.caller.set_min_qual(min_qual);
    }

    /// Drop calls deeper than `max_depth_ratio` times the median depth. 0 keeps all. Call before reading the VCF.
    pub fn set_max_depth_ratio(&mut self, max_depth_ratio: f64) {
        self.caller.set_max_depth_ratio(max_depth_ratio);
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

    /// The mean depth in windows of `window` bases as BED text with four columns (like `mosdepth` regions).
    pub fn depth_bed(&self, window: usize) -> String {
        self.caller.depth_bed(window)
    }

    /// Reads mapped so far.
    pub fn read_count(&self) -> u64 {
        self.reads
    }

    /// The weak signals so far (calls at relaxed thresholds) as VCF text. See `Caller::evidence_calls`.
    pub fn evidence_vcf(&self) -> Result<String, String> {
        let mut out = Vec::new();
        write_vcf(&mut out, &self.contigs, &self.sample, &self.caller.evidence_calls(), self.config.large_indel_min)
            .map_err(|e| e.to_string())?;
        String::from_utf8(out).map_err(|e| e.to_string())
    }

    /// The calls so far as VCF text.
    pub fn vcf(&self) -> Result<String, String> {
        let mut out = Vec::new();
        write_vcf(&mut out, &self.contigs, &self.sample, &self.caller.calls(), self.config.large_indel_min)
            .map_err(|e| e.to_string())?;
        String::from_utf8(out).map_err(|e| e.to_string())
    }
}

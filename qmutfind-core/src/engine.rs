//! The whole pipeline over bytes: a reference FASTA and FASTQ chunks go in, VCF text comes out.
//! No file access and no threads. The browser build and the CLI both drive this type.

use std::collections::VecDeque;

use crate::{
    choose_pair, read_fasta, write_vcf, Aligner, Caller, CallerConfig, FastqChunker, ReadKind, RammapAligner, DEFAULT_MAX_INSERT,
};

type Read = (String, Vec<u8>);

/// The read name without the part after the first space and without a trailing `/1` or `/2`.
fn pair_name(name: &str) -> &str {
    let first = name.split_whitespace().next().unwrap_or("");
    first.strip_suffix("/1").or_else(|| first.strip_suffix("/2")).unwrap_or(first)
}

pub struct Engine {
    aligner: RammapAligner,
    caller: Caller,
    chunker: FastqChunker,
    chunker2: FastqChunker,
    /// Reads of each mate file that wait for their mate.
    waiting: [VecDeque<Read>; 2],
    max_insert: usize,
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
            chunker2: FastqChunker::new(),
            waiting: [VecDeque::new(), VecDeque::new()],
            max_insert: DEFAULT_MAX_INSERT,
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

    /// The largest insert (start of the forward mate to end of the reverse mate) of a proper pair. Call before adding reads.
    pub fn set_max_insert(&mut self, max_insert: usize) {
        self.max_insert = max_insert;
    }

    /// Add a chunk of the first mate file of a paired run. `add_mate2` takes the second file. The chunks of the two files
    /// can have any size and can come in any order; reads are mapped as pairs when both mates are present.
    pub fn add_mate1(&mut self, chunk: &[u8]) -> Result<(), String> {
        let records = self.chunker.push(chunk).map_err(|e| format!("reads (mate 1): {e}"))?;
        self.waiting[0].extend(records);
        self.map_pairs()
    }

    /// Add a chunk of the second mate file. See `add_mate1`.
    pub fn add_mate2(&mut self, chunk: &[u8]) -> Result<(), String> {
        let records = self.chunker2.push(chunk).map_err(|e| format!("reads (mate 2): {e}"))?;
        self.waiting[1].extend(records);
        self.map_pairs()
    }

    /// Mark the end of both mate files. A truncated record, or a file with more reads than the other, is an error.
    pub fn end_mates(&mut self) -> Result<(), String> {
        let r1 = self.chunker.finish().map_err(|e| format!("reads (mate 1): {e}"))?;
        let r2 = self.chunker2.finish().map_err(|e| format!("reads (mate 2): {e}"))?;
        self.waiting[0].extend(r1);
        self.waiting[1].extend(r2);
        self.map_pairs()?;
        if !self.waiting[0].is_empty() || !self.waiting[1].is_empty() {
            let (a, b) = (self.waiting[0].len(), self.waiting[1].len());
            self.waiting = [VecDeque::new(), VecDeque::new()];
            return Err(format!("the two read files have different numbers of reads: {a} reads of mate 1 and {b} reads of mate 2 have no mate"));
        }
        Ok(())
    }

    fn map_pairs(&mut self) -> Result<(), String> {
        while !self.waiting[0].is_empty() && !self.waiting[1].is_empty() {
            let (n1, s1) = self.waiting[0].pop_front().expect("checked non-empty");
            let (n2, s2) = self.waiting[1].pop_front().expect("checked non-empty");
            if pair_name(&n1) != pair_name(&n2) {
                return Err(format!("the read files are not in the same order: mate 1 read {n1:?} is paired with mate 2 read {n2:?}"));
            }
            let (a1, a2) = choose_pair(self.aligner.align_candidates(&n1, &s1), self.aligner.align_candidates(&n2, &s2), self.max_insert);
            for alignment in &a1 {
                self.caller.add(alignment, &s1);
            }
            for alignment in &a2 {
                self.caller.add(alignment, &s2);
            }
            self.reads += 2;
        }
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

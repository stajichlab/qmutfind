//! Haploid pileup caller. Counts bases per reference position and calls sites where
//! one alternate allele has enough reads.

use std::collections::HashMap;

use crate::{Alignment, CigarOp};

/// Calling thresholds. Provisional defaults are `UNMEASURED` (see PLAN.md section 5).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CallerConfig {
    pub min_depth: u32,
    /// Minimum alternate fraction for a SNP.
    pub min_snp_fraction: f64,
    /// Minimum alternate fraction for an indel.
    pub min_indel_fraction: f64,
    /// Insertions of this length or more are reported as symbolic `<INS>` calls.
    pub large_indel_min: usize,
}

/// One called variant. `pos` is 1-based.
#[derive(Debug, Clone, PartialEq)]
pub struct Call {
    pub contig: String,
    pub pos: u64,
    pub reference: String,
    pub alt: String,
    pub depth: u32,
    pub alt_count: u32,
    /// Set for a symbolic call (`alt` is `<INS>`): the inserted length.
    pub svlen: Option<usize>,
    /// Phred-scaled probability that this many alternate reads come from sequencing error alone (a binomial tail
    /// at the error rate estimated from the whole sample). It does not cover mapping errors: a call in a repeat
    /// can have a high QUAL. `None` when the caller did not compute it.
    pub qual: Option<f64>,
    /// The depth at the call divided by the median depth of the sample. Near 1 is typical. Well above 1 suggests
    /// reads from more than one copy of a repeat. `None` when not computed.
    pub depth_ratio: Option<f64>,
}

/// The mean depth of one window of a contig. `start` is 0-based and `end` is exclusive.
#[derive(Debug, Clone, PartialEq)]
pub struct DepthWindow {
    pub contig: String,
    pub start: usize,
    pub end: usize,
    pub mean: f64,
}

/// Natural log of the gamma function (Lanczos approximation, g = 7).
fn ln_gamma(x: f64) -> f64 {
    const C: [f64; 9] = [
        0.999_999_999_999_809_9, 676.520_368_121_885_1, -1_259.139_216_722_402_8, 771.323_428_777_653_1,
        -176.615_029_162_140_6, 12.507_343_278_686_905, -0.138_571_095_265_720_12, 9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + 7.5;
    for (i, c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// Phred-scaled probability P(X >= k) for X ~ Binomial(n, e), at most `QUAL_MAX`.
fn phred_binomial_tail(n: u32, k: u32, e: f64) -> f64 {
    let (n, k) = (n as f64, k as f64);
    let (ln_e, ln_q) = (e.ln(), (1.0 - e).ln());
    let ln_choose = |i: f64| ln_gamma(n + 1.0) - ln_gamma(i + 1.0) - ln_gamma(n - i + 1.0);
    let mut terms: Vec<f64> = Vec::new();
    let mut i = k;
    while i <= n {
        terms.push(ln_choose(i) + i * ln_e + (n - i) * ln_q);
        i += 1.0;
    }
    let max = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let ln_tail = max + terms.iter().map(|t| (t - max).exp()).sum::<f64>().ln();
    (-10.0 * ln_tail / std::f64::consts::LN_10).clamp(0.0, QUAL_MAX)
}

/// The largest QUAL written. A binomial tail below this is far beyond what the error model can support.
pub const QUAL_MAX: f64 = 999.0;

/// What `calls` learns from the whole sample: the median depth and the error rates that QUAL uses.
struct SampleStats {
    median_depth: f64,
    snp_error: f64,
    indel_error: f64,
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
    min_mapq: u8,
    min_indel_alt_reads: u32,
    /// Calls with a lower QUAL are dropped. 0 keeps all.
    min_qual: f64,
    /// Calls with a larger depth ratio are dropped. 0 means no limit.
    max_depth_ratio: f64,
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

impl CallerConfig {
    /// Default minimum number of reads for an indel. For Nanopore reads, 1 bp deletions in homopolymers
    /// pass the fraction threshold at low depth. A floor of 6 reads cut the false small indels from 42 to 16
    /// (seed 1, 27.7x) and from 234 to 72 (seed 2, 15.8x) and cost 0 and 6 true indels (PLAN.md section 11).
    /// Short reads stay at 0: a floor of 6 would remove every call at 7x coverage.
    pub fn default_min_indel_alt_reads(kind: crate::ReadKind) -> u32 {
        match kind {
            crate::ReadKind::Short => 0,
            crate::ReadKind::Long => 6,
        }
    }

    /// Default minimum mapping quality by read kind. For short reads, 20 cut the EVOL20-only calls from
    /// 61 to 9 and kept both candidate variants, at no cost on the simulated reads (PLAN.md section 11).
    /// Long reads stay at 0 because nothing was measured for them.
    pub fn default_min_mapq(kind: crate::ReadKind) -> u8 {
        match kind {
            crate::ReadKind::Short => 20,
            crate::ReadKind::Long => 0,
        }
    }

    /// Default thresholds by read kind. They are set from the benchmark and the EVOL20 real data
    /// (PLAN.md section 5): a haploid SNP should be near fraction 1, and long reads have lower
    /// fractions because of their errors. Indels use a lower threshold than SNPs.
    pub fn defaults_for(kind: crate::ReadKind) -> Self {
        CallerConfig {
            min_depth: 5,
            min_snp_fraction: match kind {
                crate::ReadKind::Short => 0.8,
                crate::ReadKind::Long => 0.7,
            },
            min_indel_fraction: 0.5,
            large_indel_min: 50,
        }
    }
}

impl Caller {
    pub fn new(reference: Vec<(String, Vec<u8>)>, config: CallerConfig) -> Self {
        let counts = reference.iter().map(|(_, s)| vec![[0u16; 4]; s.len()]).collect();
        Caller { reference, counts, deletions: HashMap::new(), insertions: HashMap::new(), config, min_mapq: 0, min_indel_alt_reads: 0, min_qual: 0.0, max_depth_ratio: 0.0 }
    }

    /// Ignore alignments with a mapping quality below `min_mapq`. The default is 0, which keeps all.
    pub fn set_min_mapq(&mut self, min_mapq: u8) {
        self.min_mapq = min_mapq;
    }

    /// The minimum mapping quality in use.
    pub fn min_mapq(&self) -> u8 {
        self.min_mapq
    }

    /// Call an indel only when at least this many reads support it. The default is 0, which keeps all.
    pub fn set_min_indel_alt_reads(&mut self, min_alt_reads: u32) {
        self.min_indel_alt_reads = min_alt_reads;
    }

    /// The minimum number of reads for an indel in use.
    pub fn min_indel_alt_reads(&self) -> u32 {
        self.min_indel_alt_reads
    }

    /// Drop calls with a QUAL below `min_qual`. The default is 0, which keeps all.
    pub fn set_min_qual(&mut self, min_qual: f64) {
        self.min_qual = min_qual;
    }

    /// The minimum QUAL in use.
    pub fn min_qual(&self) -> f64 {
        self.min_qual
    }

    /// Drop calls whose depth is more than `max_depth_ratio` times the median depth. 0 (the default) keeps all.
    pub fn set_max_depth_ratio(&mut self, max_depth_ratio: f64) {
        self.max_depth_ratio = max_depth_ratio;
    }

    /// The maximum depth ratio in use.
    pub fn max_depth_ratio(&self) -> f64 {
        self.max_depth_ratio
    }

    /// Replace the thresholds. Call before `calls`.
    pub fn set_config(&mut self, config: CallerConfig) {
        self.config = config;
    }

    /// Add one alignment. `read` is the sequence that was given to the aligner.
    pub fn add(&mut self, alignment: &Alignment, read: &[u8]) {
        if alignment.mapq() < self.min_mapq {
            return;
        }
        let aligned = alignment.aligned_read(read);
        self.add_parts(&alignment.contig, alignment.ref_start, alignment.cigar_ops(), &aligned);
    }

    /// Add one record read from a SAM file.
    pub fn add_sam(&mut self, record: &crate::SamRecord) {
        if record.mapq < self.min_mapq {
            return;
        }
        self.add_parts(&record.contig, record.ref_start, &record.cigar_ops, &record.seq);
    }

    /// Walk the CIGAR operations and count the bases. `aligned` holds the aligned bases in
    /// reference orientation, without clipped bases.
    fn add_parts(&mut self, contig: &str, ref_start: usize, ops: &[CigarOp], aligned: &[u8]) {
        let Some(ci) = self.reference.iter().position(|(n, _)| n == contig) else {
            return;
        };
        let (mut qi, mut ri) = (0usize, ref_start);
        for op in ops {
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

    /// Mean depth in consecutive windows of `window` bases along each contig, uncovered bases counted as zero
    /// (the last window of a contig may be shorter). Like `mosdepth` with `--by window`.
    pub fn depth_windows(&self, window: usize) -> Vec<DepthWindow> {
        let window = window.max(1);
        let mut out = Vec::new();
        for (ci, (name, seq)) in self.reference.iter().enumerate() {
            let mut start = 0;
            while start < seq.len() {
                let end = (start + window).min(seq.len());
                let sum: u64 = self.counts[ci][start..end].iter().map(|c| c.iter().map(|&x| x as u64).sum::<u64>()).sum();
                out.push(DepthWindow { contig: name.clone(), start, end, mean: sum as f64 / (end - start) as f64 });
                start = end;
            }
        }
        out
    }

    /// The read depth at each position of `contig` from 1-based `start` to `end` (inclusive), or None when there is no such contig.
    pub(crate) fn depth_range(&self, contig: &str, start: u64, end: u64) -> Option<Vec<u32>> {
        let ci = self.reference.iter().position(|(n, _)| n == contig)?;
        let len = self.counts[ci].len() as u64;
        let (lo, hi) = (start.max(1), end.min(len));
        Some((lo..=hi).map(|p| self.counts[ci][(p - 1) as usize].iter().map(|&x| x as u32).sum()).collect())
    }

    /// A histogram of the depth (0 to 2,000, higher values counted as 2,000) over the positions of every contig not in `exclude` that have
    /// at least one read.
    pub(crate) fn depth_histogram(&self, exclude: &[String]) -> Vec<u64> {
        let mut histogram = vec![0u64; 2_001];
        for (ci, (name, _)) in self.reference.iter().enumerate() {
            if exclude.contains(name) {
                continue;
            }
            for c in &self.counts[ci] {
                let depth = c.iter().map(|&x| x as usize).sum::<usize>();
                if depth > 0 {
                    histogram[depth.min(2_000)] += 1;
                }
            }
        }
        histogram
    }

    /// `depth_windows` as BED text with four columns (contig, start, end, mean depth), the format of the
    /// `mosdepth` regions file.
    pub fn depth_bed(&self, window: usize) -> String {
        self.depth_windows(window).iter().map(|w| format!("{}\t{}\t{}\t{:.2}\n", w.contig, w.start, w.end, w.mean)).collect()
    }

    /// The median depth over covered positions, the rate of non-reference bases (SNP error) and the rate of indel
    /// events per aligned base (indel error). The rates have a small floor so that QUAL stays finite.
    fn sample_stats(&self) -> SampleStats {
        let mut histogram = vec![0u64; 2_001];
        let (mut covered, mut total_depth, mut non_reference) = (0u64, 0u64, 0u64);
        for (ci, (_, seq)) in self.reference.iter().enumerate() {
            for (pos, c) in self.counts[ci].iter().enumerate() {
                let depth: u64 = c.iter().map(|&x| x as u64).sum();
                if depth == 0 {
                    continue;
                }
                covered += 1;
                total_depth += depth;
                histogram[(depth as usize).min(2_000)] += 1;
                if let Some(r) = base_index(seq[pos]) {
                    non_reference += depth - c[r] as u64;
                }
            }
        }
        let mut seen = 0u64;
        let mut median = 0usize;
        for (d, &n) in histogram.iter().enumerate() {
            seen += n;
            if covered > 0 && seen * 2 >= covered {
                median = d;
                break;
            }
        }
        let indel_events: u64 = self.deletions.values().map(|&n| n as u64).sum::<u64>() + self.insertions.values().map(|e| e.count as u64).sum::<u64>();
        let rate = |events: u64, floor: f64| if total_depth == 0 { floor } else { (events as f64 / total_depth as f64).max(floor) };
        SampleStats { median_depth: median as f64, snp_error: rate(non_reference, 1e-6), indel_error: rate(indel_events, 1e-7) }
    }

    /// All calls, in reference order.
    pub fn calls(&self) -> Vec<Call> {
        self.calls_with(&self.config, self.min_indel_alt_reads, true)
    }

    /// Calls at relaxed thresholds (depth 2, fraction 0.1, two reads for an indel, no QUAL or depth-ratio filter): the
    /// weak signals below the calling thresholds. `combine` uses them to tell "some reads carry it" from "none do".
    /// It is a superset of `calls()` for SNPs; an indel needs two supporting reads here, as in the strict set only
    /// when `min_indel_alt_reads` is lower.
    pub fn evidence_calls(&self) -> Vec<Call> {
        let relaxed = CallerConfig { min_depth: 2, min_snp_fraction: 0.1, min_indel_fraction: 0.1, large_indel_min: self.config.large_indel_min };
        self.calls_with(&relaxed, 2, false)
    }

    fn calls_with(&self, config: &CallerConfig, min_indel_alt_reads: u32, filters: bool) -> Vec<Call> {
        let stats = self.sample_stats();
        let depth_ratio = |depth: u32| if stats.median_depth > 0.0 { Some(depth as f64 / stats.median_depth) } else { None };
        let snp_qual = |depth: u32, alt: u32| Some(phred_binomial_tail(depth, alt, stats.snp_error));
        let indel_qual = |depth: u32, alt: u32| Some(phred_binomial_tail(depth, alt, stats.indel_error));
        let mut out: Vec<(usize, Call)> = Vec::new();
        for (ci, (name, seq)) in self.reference.iter().enumerate() {
            for (pos, c) in self.counts[ci].iter().enumerate() {
                let depth: u32 = c.iter().map(|&x| x as u32).sum();
                if depth < config.min_depth {
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
                if alt_count > 0 && alt_count as f64 / depth as f64 >= config.min_snp_fraction {
                    out.push((ci, Call {
                        contig: name.clone(),
                        pos: pos as u64 + 1,
                        reference: (seq[pos] as char).to_string(),
                        alt: (b"ACGT"[best] as char).to_string(),
                        depth,
                        alt_count,
                        svlen: None,
                        qual: snp_qual(depth, alt_count),
                        depth_ratio: depth_ratio(depth),
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
            if depth >= config.min_depth
                && alt_count as f64 / depth as f64 >= config.min_indel_fraction
                && alt_count >= min_indel_alt_reads
            {
                out.push((ci, Call {
                    contig: self.reference[ci].0.clone(),
                    pos: anchor as u64 + 1,
                    reference: String::from_utf8_lossy(&seq[anchor..anchor + 1 + len]).into_owned(),
                    alt: (seq[anchor] as char).to_string(),
                    depth,
                    alt_count,
                    svlen: None,
                    qual: indel_qual(depth, alt_count),
                    depth_ratio: depth_ratio(depth),
                }));
            }
        }
        let mut large: Vec<Vec<(usize, usize, u32)>> = vec![Vec::new(); self.reference.len()];
        for (&(ci, anchor, len), ev) in &self.insertions {
            if len >= config.large_indel_min {
                large[ci].push((anchor, len, ev.count));
                continue;
            }
            let alt_count = ev.count;
            let seq = &self.reference[ci].1;
            if base_index(seq[anchor]).is_none() {
                continue; // anchor base is ambiguous
            }
            let depth: u32 = self.counts[ci][anchor].iter().map(|&x| x as u32).sum();
            if depth >= config.min_depth
                && alt_count as f64 / depth as f64 >= config.min_indel_fraction
                && alt_count >= min_indel_alt_reads
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
                    svlen: None,
                    qual: indel_qual(depth, alt_count),
                    depth_ratio: depth_ratio(depth),
                }));
            }
        }
        for (ci, events) in large.into_iter().enumerate() {
            let seq = &self.reference[ci].1;
            for (anchor, svlen, support) in cluster_large_insertions(events) {
                if base_index(seq[anchor]).is_none() {
                    continue; // anchor base is ambiguous
                }
                let depth: u32 = self.counts[ci][anchor].iter().map(|&x| x as u32).sum();
                if depth >= config.min_depth
                    && support as f64 / depth as f64 >= config.min_indel_fraction
                    && support >= min_indel_alt_reads
                {
                    out.push((ci, Call {
                        contig: self.reference[ci].0.clone(),
                        pos: anchor as u64 + 1,
                        reference: (seq[anchor] as char).to_string(),
                        alt: "<INS>".to_string(),
                        depth,
                        alt_count: support,
                        svlen: Some(svlen),
                        qual: indel_qual(depth, support),
                        depth_ratio: depth_ratio(depth),
                    }));
                }
            }
        }
        out.sort_by_key(|(ci, c)| (*ci, c.pos));
        let kept: Vec<Call> = out
            .into_iter()
            .map(|(_, c)| c)
            .filter(|c| !filters || c.qual.map_or(true, |q| q >= self.min_qual))
            .filter(|c| !filters || self.max_depth_ratio <= 0.0 || c.depth_ratio.map_or(true, |r| r <= self.max_depth_ratio))
            .collect();
        self.left_aligned(kept)
    }

    /// Write every explicit indel in its leftmost form, so that the same event has the same REF, ALT and POS in every
    /// sample. Reads place an indel in a repeat at different bases; the evidence is counted where the reads put it.
    /// Calls that become equal keep the one with more supporting reads (the counts are not added). The result is
    /// sorted by contig order and position.
    fn left_aligned(&self, calls: Vec<Call>) -> Vec<Call> {
        let rank: HashMap<&str, usize> = self.reference.iter().enumerate().map(|(i, (name, _))| (name.as_str(), i)).collect();
        let mut aligned = crate::normalize::left_align_calls(calls, &self.reference);
        aligned.sort_by_key(|c| (rank.get(c.contig.as_str()).copied().unwrap_or(usize::MAX), c.pos));
        let mut out: Vec<Call> = Vec::with_capacity(aligned.len());
        for c in aligned {
            match out.last_mut() {
                Some(last) if last.contig == c.contig && last.pos == c.pos && last.reference == c.reference && last.alt == c.alt => {
                    if c.alt_count > last.alt_count {
                        *last = c;
                    }
                }
                _ => out.push(c),
            }
        }
        out
    }
}

/// Position tolerance, in bases, for events that belong to one large insertion.
const ANCHOR_TOLERANCE: usize = 20;

/// Length tolerance for events that belong to one large insertion: 10% of the length, at least 5.
fn length_tolerance(len: usize) -> usize {
    (len / 10).max(5)
}

/// Group large insertion events (anchor, length, reads) by position and length.
/// Returns one (anchor, median length, supporting reads) per group. The anchor is the
/// most supported one in the group.
fn cluster_large_insertions(mut events: Vec<(usize, usize, u32)>) -> Vec<(usize, usize, u32)> {
    events.sort();
    let mut out = Vec::new();
    let mut i = 0;
    while i < events.len() {
        let mut j = i + 1;
        while j < events.len() && events[j].0 - events[j - 1].0 <= ANCHOR_TOLERANCE {
            j += 1;
        }
        let mut near: Vec<(usize, usize, u32)> = events[i..j].to_vec();
        near.sort_by_key(|e| e.1);
        let mut k = 0;
        while k < near.len() {
            let mut m = k + 1;
            while m < near.len() && near[m].1 - near[m - 1].1 <= length_tolerance(near[m - 1].1) {
                m += 1;
            }
            let group = &near[k..m];
            let support: u32 = group.iter().map(|e| e.2).sum();
            let anchor = group
                .iter()
                .max_by_key(|e| (e.2, std::cmp::Reverse(e.0)))
                .map(|e| e.0)
                .unwrap();
            let mut seen = 0;
            let mut median = group[0].1;
            for e in group {
                seen += e.2;
                if seen * 2 >= support {
                    median = e.1;
                    break;
                }
            }
            out.push((anchor, median, support));
            k = m;
        }
        i = j;
    }
    out
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

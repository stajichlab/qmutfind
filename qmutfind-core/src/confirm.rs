//! Confirm an introduced allele: a marker cassette (for example `hph`) that replaces a gene or is inserted somewhere. The reference is the
//! unmodified genome plus the cassette as one more contig; the reads come from the transformant. From the alignments the module answers
//! four questions: where does the cassette join the genome (junctions of split reads), is that the intended place, is the cassette intact
//! (coverage, variants, internal junctions), and how many copies are there (cassette depth against genome depth).
//!
//! A junction is a read whose two parts align to different places, next to each other on the read. Each end of a junction is the last
//! reference base of one part (`pos`, 1-based) and the side of that base on which the aligned part lies (`Side`). Ends are written so that a
//! read from either strand gives the same junction. The checks give Pass, Fail or Inconclusive (not enough coverage to tell); they never turn
//! a lack of reads into a pass. The thresholds are in `ConfirmOptions` and are not tuned on real data.

use crate::caller::Call;
use crate::{Alignment, Engine, ReadKind};
use std::collections::HashMap;

/// Which side of a breakpoint the aligned part of a read lies on, in reference coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Side {
    Left,
    Right,
}

/// One end of a junction: the reference base `pos` (1-based) next to the breakpoint, and the side where the aligned bases are.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct End {
    pub contig: String,
    pub pos: u64,
    pub keeps: Side,
}

/// A group of split reads with the same two ends (within a small tolerance).
#[derive(Debug, Clone, PartialEq)]
pub struct Junction {
    pub a: End,
    pub b: End,
    pub reads: u32,
    /// Bases that both parts of the read could explain (positive: the two ends share this many bases) or bases that neither explains (negative: inserted).
    pub overlap: i32,
}

/// How the cassette lies in the genome: `Forward` when its start is joined to the left flank.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Forward,
    Reverse,
}

/// The place where a cassette is joined to the genome: the two junctions of one insertion (either can be missing).
#[derive(Debug, Clone, PartialEq)]
pub struct Site {
    pub contig: String,
    pub construct: String,
    /// The last genome base before the cassette (the left flank), when that junction was found.
    pub left: Option<u64>,
    /// The first genome base after the cassette (the right flank).
    pub right: Option<u64>,
    pub reads_left: u32,
    pub reads_right: u32,
    pub orientation: Option<Orientation>,
    /// Genome bases between the two flanks, which the insertion replaces (0 for a plain insertion).
    pub replaced_bases: Option<i64>,
    /// Where on the cassette each junction lies (1-based).
    pub construct_at_left: Option<u64>,
    pub construct_at_right: Option<u64>,
    /// Bases that the genome and the cassette share at each junction (0 when none; negative when bases are inserted between them). With shared
    /// bases the breakpoint could be that many bases further out; `left` and `right` are the bases that are not in doubt.
    pub homology_left: i32,
    pub homology_right: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckStatus {
    Pass,
    Fail,
    /// Not enough coverage to tell.
    Inconclusive,
    NotApplicable,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Check {
    pub name: &'static str,
    pub status: CheckStatus,
    pub detail: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// No check failed and none was inconclusive.
    Confirmed,
    /// At least one check failed.
    NotConfirmed,
    /// No check failed, but at least one could not be decided.
    Inconclusive,
}

/// What the reads say about one construct contig.
#[derive(Debug, Clone, PartialEq)]
pub struct CassetteSummary {
    pub name: String,
    pub length: u64,
    pub median_depth: f64,
    /// Median depth of the cassette divided by the median depth of the genome: about the number of copies per genome.
    pub copy_ratio: f64,
    /// The read depth at each position of the cassette (the first value is position 1).
    pub depth: Vec<u32>,
    /// Stretches (1-based, inclusive) with too little coverage.
    pub gaps: Vec<(u64, u64)>,
    /// Calls on the cassette: differences from the sequence that was given.
    pub variants: Vec<Call>,
}

/// The region that the cassette should replace (1-based, inclusive).
#[derive(Debug, Clone, PartialEq)]
pub struct Target {
    pub contig: String,
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Clone)]
pub struct ConfirmOptions {
    pub target: Option<Target>,
    /// A junction needs this many reads.
    pub min_junction_reads: u32,
    /// Each part of a split read needs this many bases and this mapping quality.
    pub min_part_len: usize,
    pub min_part_mapq: u8,
    /// Junction ends within this many bases of each other are one junction.
    pub slop: u64,
    /// A junction is at the target when each flank is within this many bases of the target edge.
    pub target_tolerance: u64,
    /// Below this genome depth, a check that rests on missing reads is inconclusive.
    pub min_depth: f64,
    /// Two junctions of one insertion are at most this far apart in the genome.
    pub max_site_span: u64,
    /// A cassette base has a gap when its depth is below `max(3, gap_fraction * genome depth)`.
    pub gap_fraction: f64,
}

impl Default for ConfirmOptions {
    fn default() -> Self {
        ConfirmOptions { target: None, min_junction_reads: 3, min_part_len: 25, min_part_mapq: 1, slop: 3, target_tolerance: 20, min_depth: 10.0, max_site_span: 200_000, gap_fraction: 0.2 }
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct JunctionOptions {
    pub min_part_len: usize,
    pub min_part_mapq: u8,
}

/// Collects the junctions of split reads while the reads are mapped.
pub(crate) struct JunctionCollector {
    order: HashMap<String, usize>,
    options: JunctionOptions,
    observed: Vec<(End, End, i32)>,
}

/// Two parts are neighbours on the read when they overlap by at most this many bases, or are separated by at most this many.
const MAX_OVERLAP: i64 = 50;
const MAX_INSERT: i64 = 50;
/// Two parts on one contig, one strand and this close are a gapped alignment, not a junction.
const MIN_SPLIT_DISTANCE: u64 = 50;

impl JunctionCollector {
    pub(crate) fn new(order: Vec<String>, options: JunctionOptions) -> Self {
        JunctionCollector { order: order.into_iter().enumerate().map(|(i, n)| (n, i)).collect(), options, observed: Vec::new() }
    }

    /// The end of the part where the read leaves it (the later end on the read).
    fn exit_end(a: &Alignment) -> End {
        if a.reverse {
            End { contig: a.contig.clone(), pos: a.ref_start as u64 + 1, keeps: Side::Right }
        } else {
            End { contig: a.contig.clone(), pos: a.ref_end() as u64, keeps: Side::Left }
        }
    }

    /// The end of the part where the read enters it (the earlier end on the read).
    fn entry_end(a: &Alignment) -> End {
        if a.reverse {
            End { contig: a.contig.clone(), pos: a.ref_end() as u64, keeps: Side::Left }
        } else {
            End { contig: a.contig.clone(), pos: a.ref_start as u64 + 1, keeps: Side::Right }
        }
    }

    pub(crate) fn add(&mut self, alignments: &[Alignment]) {
        if alignments.len() < 2 {
            return;
        }
        let mut parts: Vec<&Alignment> = alignments
            .iter()
            .filter(|a| {
                let (s, e) = a.query_range();
                e.saturating_sub(s) >= self.options.min_part_len && a.mapq() >= self.options.min_part_mapq
            })
            .collect();
        if parts.len() < 2 {
            return;
        }
        parts.sort_by_key(|a| a.query_range().0);
        for pair in parts.windows(2) {
            let (x, y) = (pair[0], pair[1]);
            let overlap = x.query_range().1 as i64 - y.query_range().0 as i64;
            if overlap > MAX_OVERLAP || -overlap > MAX_INSERT {
                continue;
            }
            if x.contig == y.contig && x.reverse == y.reverse && (x.ref_start as u64).abs_diff(y.ref_start as u64) < MIN_SPLIT_DISTANCE {
                continue;
            }
            let (mut a, mut b) = (Self::exit_end(x), Self::entry_end(y));
            // The same junction from a read of the other strand has its ends in the other order; sort them.
            let key = |e: &End| (self.order.get(&e.contig).copied().unwrap_or(usize::MAX), e.pos, e.keeps);
            if key(&a) > key(&b) {
                std::mem::swap(&mut a, &mut b);
            }
            self.observed.push((a, b, overlap as i32));
        }
    }

    /// The junctions with at least `min_reads` reads. Observations whose ends are within `slop` bases of each other are one junction.
    pub(crate) fn junctions(&self, slop: u64, min_reads: u32) -> Vec<Junction> {
        let mut groups: HashMap<(String, Side, String, Side), Vec<&(End, End, i32)>> = HashMap::new();
        for o in &self.observed {
            groups.entry((o.0.contig.clone(), o.0.keeps, o.1.contig.clone(), o.1.keeps)).or_default().push(o);
        }
        let mut out = Vec::new();
        for (_, mut members) in groups {
            members.sort_by_key(|o| (o.0.pos, o.1.pos));
            let mut clusters: Vec<Vec<&(End, End, i32)>> = Vec::new();
            for o in members {
                match clusters.iter_mut().find(|c| c[0].0.pos.abs_diff(o.0.pos) <= slop && c[0].1.pos.abs_diff(o.1.pos) <= slop) {
                    Some(c) => c.push(o),
                    None => clusters.push(vec![o]),
                }
            }
            for c in clusters {
                if (c.len() as u32) < min_reads {
                    continue;
                }
                // The most common pair of positions is the junction; the median overlap says whether bases are shared or inserted.
                let mut counts: HashMap<(u64, u64), u32> = HashMap::new();
                for o in &c {
                    *counts.entry((o.0.pos, o.1.pos)).or_insert(0) += 1;
                }
                let ((pa, pb), _) = counts.into_iter().max_by_key(|&((a, b), n)| (n, std::cmp::Reverse((a, b)))).expect("a cluster has members");
                let mut overlaps: Vec<i32> = c.iter().map(|o| o.2).collect();
                overlaps.sort();
                out.push(Junction {
                    a: End { pos: pa, ..c[0].0.clone() },
                    b: End { pos: pb, ..c[0].1.clone() },
                    reads: c.len() as u32,
                    overlap: overlaps[overlaps.len() / 2],
                });
            }
        }
        out.sort_by(|x, y| y.reads.cmp(&x.reads).then_with(|| x.a.cmp(&y.a)).then_with(|| x.b.cmp(&y.b)));
        out
    }
}

/// Length of the words used to find construct sequence that is also in the genome.
const SHARED_WORD: usize = 24;

/// A warning for each construct that shares sequence with the genome (24-base words found in both, on either strand): the usual cause is homology arms
/// in the construct. Reads from such a stretch map to the genome and to the construct equally well, so coverage, junctions and copy number of the
/// cassette cannot be trusted; the construct should hold the cassette only.
fn shared_with_genome(reference_fasta: &[u8], constructs: &[&str]) -> Result<Vec<String>, String> {
    use std::collections::HashSet;
    let contigs = crate::read_fasta(reference_fasta).map_err(|e| format!("reference FASTA: {e}"))?;
    let code = |b: u8| match b.to_ascii_uppercase() { b'A' => Some(0u64), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None };
    let mask = (1u64 << (2 * SHARED_WORD)) - 1;
    // Every word of the constructs, forward and reverse complement, with the construct and position of the forward word.
    let mut words: std::collections::HashMap<u64, Vec<(usize, usize)>> = std::collections::HashMap::new();
    let mut construct_index: Vec<(&str, usize)> = Vec::new();
    for (ci, (name, seq)) in contigs.iter().enumerate() {
        if !constructs.contains(&name.as_str()) { continue }
        construct_index.push((name.as_str(), ci));
        let (mut fwd, mut rev, mut valid) = (0u64, 0u64, 0usize);
        for (i, &b) in seq.iter().enumerate() {
            match code(b) {
                Some(c) => {
                    fwd = ((fwd << 2) | c) & mask;
                    rev = (rev >> 2) | ((3 - c) << (2 * (SHARED_WORD - 1)));
                    valid += 1;
                }
                None => valid = 0,
            }
            if valid >= SHARED_WORD {
                let start = i + 1 - SHARED_WORD;
                words.entry(fwd).or_default().push((ci, start));
                words.entry(rev).or_default().push((ci, start));
            }
        }
    }
    let mut shared: HashSet<(usize, usize)> = HashSet::new(); // (construct, word start) found in the genome
    let own: HashSet<&str> = constructs.iter().copied().collect();
    for (name, seq) in &contigs {
        if own.contains(name.as_str()) { continue }
        let (mut fwd, mut valid) = (0u64, 0usize);
        for &b in seq {
            match code(b) {
                Some(c) => { fwd = ((fwd << 2) | c) & mask; valid += 1 }
                None => valid = 0,
            }
            if valid >= SHARED_WORD {
                if let Some(hits) = words.get(&fwd) {
                    shared.extend(hits.iter().copied());
                }
            }
        }
    }
    let mut warnings = Vec::new();
    for (name, ci) in construct_index {
        let mut starts: Vec<usize> = shared.iter().filter(|(c, _)| *c == ci).map(|&(_, s)| s).collect();
        starts.sort_unstable();
        // Merge the words into stretches of construct positions (1-based, inclusive).
        let mut stretches: Vec<(usize, usize)> = Vec::new();
        for s in starts {
            match stretches.last_mut() {
                Some(last) if s <= last.1 => last.1 = s + SHARED_WORD,
                _ => stretches.push((s + 1, s + SHARED_WORD)),
            }
        }
        let total: usize = stretches.iter().map(|(a, b)| b - a + 1).sum();
        if total >= 50 {
            let list: Vec<String> = stretches.iter().filter(|(a, b)| b - a + 1 >= 30).take(6).map(|(a, b)| format!("{a}-{b}")).collect();
            warnings.push(format!(
                "{total} bases of the construct {name} also occur in the genome (positions {}). If these are homology arms, reads from them map to the genome and to the construct equally well, so the coverage, junctions and copy number of the cassette cannot be trusted. Give the cassette sequence without the arms.",
                list.join(", ")
            ));
        }
    }
    Ok(warnings)
}

/// The region of a gene in a GFF3 file, found by its `ID` (with or without the `gene-` of NCBI files), `Name`, `gene` or `locus_tag` (or `old_locus_tag`).
/// Only `gene` and `pseudogene` lines count: the name of an mRNA is not a gene.
pub fn find_gene<R: std::io::BufRead>(gff: R, name: &str) -> Result<Target, String> {
    for line in gff.lines() {
        let line = line.map_err(|e| format!("GFF: {e}"))?;
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || !matches!(f[2], "gene" | "pseudogene") {
            continue;
        }
        let hit = f[8].split(';').filter_map(|kv| kv.trim().split_once('=')).any(|(k, v)| {
            matches!(k, "ID" | "Name" | "gene" | "locus_tag" | "old_locus_tag") && (v == name || (k == "ID" && v.strip_prefix("gene-") == Some(name)) || (k == "ID" && name.strip_prefix("gene-") == Some(v)))
        });
        if hit {
            let start: u64 = f[3].parse().map_err(|_| format!("GFF: bad start in {line}"))?;
            let end: u64 = f[4].parse().map_err(|_| format!("GFF: bad end in {line}"))?;
            return Ok(Target { contig: f[0].to_string(), start, end });
        }
    }
    Err(format!("the GFF has no gene named {name} (ID, Name, gene, locus_tag or old_locus_tag of a gene line)"))
}

/// Runs the mapping and gives the report. Reference: the genome contigs, then the construct contigs.
pub struct Confirmer {
    engine: Engine,
    constructs: Vec<String>,
    options: ConfirmOptions,
    warnings: Vec<String>,
}

impl Confirmer {
    /// `reference_fasta` holds the genome and the construct contigs (one FASTA); `constructs` names the construct contigs.
    pub fn new(reference_fasta: &[u8], constructs: &[&str], kind: ReadKind, options: ConfirmOptions) -> Result<Confirmer, String> {
        let mut engine = Engine::new(reference_fasta, kind, "confirm")?;
        for name in constructs {
            if !engine.contigs.iter().any(|(n, _)| n == name) {
                return Err(format!("the reference has no contig named {name} for the construct"));
            }
        }
        if constructs.len() >= engine.contigs.len() {
            return Err("the reference holds only construct contigs: the genome is missing".to_string());
        }
        if let Some(t) = &options.target {
            if !engine.contigs.iter().any(|(n, _)| *n == t.contig) || constructs.contains(&t.contig.as_str()) {
                return Err(format!("the target contig {} is not a genome contig of the reference", t.contig));
            }
        }
        engine.track_junctions(JunctionOptions { min_part_len: options.min_part_len, min_part_mapq: options.min_part_mapq });
        let warnings = shared_with_genome(reference_fasta, constructs)?;
        Ok(Confirmer { engine, constructs: constructs.iter().map(|s| s.to_string()).collect(), options, warnings })
    }

    /// Add a chunk of FASTQ bytes (see `Engine::add_reads`).
    pub fn add_reads(&mut self, chunk: &[u8]) -> Result<(), String> {
        self.engine.add_reads(chunk)
    }

    /// Mark the end of one reads file.
    pub fn end_file(&mut self) -> Result<(), String> {
        self.engine.end_file()
    }

    pub fn read_count(&self) -> u64 {
        self.engine.read_count()
    }

    pub fn report(&self) -> Report {
        build_report(self)
    }
}

/// The result: what was found, and the checks.
#[derive(Debug, Clone)]
pub struct Report {
    pub reads: u64,
    /// Median depth over the covered genome positions (the construct contigs are left out).
    pub genome_depth: f64,
    /// Every junction with enough reads: the two ends and the reads. Includes junctions that are not cassette junctions.
    pub junctions: Vec<Junction>,
    pub sites: Vec<Site>,
    pub cassettes: Vec<CassetteSummary>,
    pub checks: Vec<Check>,
    pub target: Option<Target>,
    /// Things that make the result less reliable (for example a construct that shares sequence with the genome).
    pub warnings: Vec<String>,
    /// The mean depth inside the target (without its edges), divided by the genome depth.
    pub target_depth_ratio: Option<f64>,
}

fn is_construct(constructs: &[String], contig: &str) -> bool {
    constructs.iter().any(|c| c == contig)
}

fn build_report(c: &Confirmer) -> Report {
    let o = &c.options;
    let collector = c.engine.junctions.as_ref().expect("Confirmer tracks junctions");
    let junctions = collector.junctions(o.slop, o.min_junction_reads);
    let histogram = c.engine.caller.depth_histogram(&c.constructs);
    let covered: u64 = histogram.iter().sum();
    let mut seen = 0u64;
    let mut genome_depth = 0.0;
    for (d, &n) in histogram.iter().enumerate() {
        seen += n;
        if covered > 0 && seen * 2 >= covered {
            genome_depth = d as f64;
            break;
        }
    }
    let enough = genome_depth >= o.min_depth;
    let low = format!("the genome coverage is only {genome_depth:.1}x (at least {:.0}x is needed to rely on missing reads)", o.min_depth);

    // Sites: pair the left and right junctions of each insertion.
    let mut lefts: Vec<(&Junction, &End, &End)> = Vec::new(); // (junction, genome end, construct end)
    let mut rights: Vec<(&Junction, &End, &End)> = Vec::new();
    let mut internal: Vec<&Junction> = Vec::new();
    for j in &junctions {
        match (is_construct(&c.constructs, &j.a.contig), is_construct(&c.constructs, &j.b.contig)) {
            (false, true) | (true, false) => {
                let (g, k) = if is_construct(&c.constructs, &j.a.contig) { (&j.b, &j.a) } else { (&j.a, &j.b) };
                if g.keeps == Side::Left { lefts.push((j, g, k)) } else { rights.push((j, g, k)) }
            }
            (true, true) => internal.push(j),
            (false, false) => {}
        }
    }
    let mut sites: Vec<Site> = Vec::new();
    let mut used_right = vec![false; rights.len()];
    let orientation_of = |g_keeps: Side, k: &End| match (g_keeps, k.keeps) {
        (Side::Left, Side::Right) | (Side::Right, Side::Left) => Orientation::Forward,
        _ => Orientation::Reverse,
    };
    for (j, g, k) in &lefts {
        // The nearest right flank after this left flank, on the same contig and the same construct, within the span.
        let best = rights
            .iter()
            .enumerate()
            .filter(|(i, (_, rg, rk))| !used_right[*i] && rg.contig == g.contig && rk.contig == k.contig && rg.pos > g.pos && rg.pos - g.pos <= o.max_site_span)
            .min_by_key(|(_, (_, rg, _))| rg.pos);
        match best {
            Some((i, (rj, rg, rk))) => {
                used_right[i] = true;
                // With bases shared by the genome and the cassette, the flank that is not in doubt is that many bases in from the aligned end.
                let (left, right) = (g.pos.saturating_sub(j.overlap.max(0) as u64), rg.pos + rj.overlap.max(0) as u64);
                sites.push(Site {
                    contig: g.contig.clone(), construct: k.contig.clone(), left: Some(left), right: Some(right), reads_left: j.reads, reads_right: rj.reads,
                    orientation: Some(orientation_of(Side::Left, k)), replaced_bases: Some(right as i64 - left as i64 - 1),
                    construct_at_left: Some(k.pos), construct_at_right: Some(rk.pos), homology_left: j.overlap, homology_right: rj.overlap,
                });
            }
            None => sites.push(Site {
                contig: g.contig.clone(), construct: k.contig.clone(), left: Some(g.pos.saturating_sub(j.overlap.max(0) as u64)), right: None, reads_left: j.reads, reads_right: 0,
                orientation: Some(orientation_of(Side::Left, k)), replaced_bases: None, construct_at_left: Some(k.pos), construct_at_right: None,
                homology_left: j.overlap, homology_right: 0,
            }),
        }
    }
    for (i, (rj, rg, rk)) in rights.iter().enumerate() {
        if !used_right[i] {
            sites.push(Site {
                contig: rg.contig.clone(), construct: rk.contig.clone(), left: None, right: Some(rg.pos + rj.overlap.max(0) as u64), reads_left: 0, reads_right: rj.reads,
                orientation: Some(orientation_of(Side::Right, rk)), replaced_bases: None, construct_at_left: None, construct_at_right: Some(rk.pos),
                homology_left: 0, homology_right: rj.overlap,
            });
        }
    }
    let target = o.target.clone();
    let at_target = |s: &Site| match (&target, s.left, s.right) {
        (Some(t), Some(l), Some(r)) => s.contig == t.contig && l.abs_diff(t.start.saturating_sub(1)) <= o.target_tolerance && r.abs_diff(t.end + 1) <= o.target_tolerance,
        _ => false,
    };
    // The main site: the one at the target, else the two-sided one with most reads, else the one with most reads.
    sites.sort_by_key(|s| {
        (!at_target(s), s.left.is_none() || s.right.is_none(), std::cmp::Reverse(s.reads_left + s.reads_right), s.left.or(s.right).unwrap_or(0))
    });

    // The cassettes.
    let gap_depth = (o.gap_fraction * genome_depth).max(3.0);
    let calls = c.engine.caller.calls();
    let mut cassettes = Vec::new();
    for name in &c.constructs {
        let length = c.engine.contigs.iter().find(|(n, _)| n == name).map(|(_, l)| *l as u64).unwrap_or(0);
        let depth = c.engine.caller.depth_range(name, 1, length).unwrap_or_default();
        let mut sorted = depth.clone();
        sorted.sort_unstable();
        let median = if sorted.is_empty() { 0.0 } else { sorted[sorted.len() / 2] as f64 };
        let trim = 5usize.min(depth.len() / 4);
        let mut gaps: Vec<(u64, u64)> = Vec::new();
        let mut run: Option<u64> = None;
        for (i, &d) in depth.iter().enumerate() {
            let pos = i as u64 + 1;
            let inside = i >= trim && i < depth.len().saturating_sub(trim);
            let gap = inside && (d as f64) < gap_depth;
            match (gap, run) {
                (true, None) => run = Some(pos),
                (false, Some(start)) => { gaps.push((start, pos - 1)); run = None }
                _ => {}
            }
        }
        if let Some(start) = run {
            gaps.push((start, depth.len() as u64 - trim as u64));
        }
        let variants: Vec<Call> = calls.iter().filter(|v| &v.contig == name).cloned().collect();
        cassettes.push(CassetteSummary { name: name.clone(), length, median_depth: median, copy_ratio: if genome_depth > 0.0 { median / genome_depth } else { 0.0 }, depth, gaps, variants });
    }
    let any_cassette_reads = cassettes.iter().any(|k| k.median_depth >= gap_depth);

    // The depth inside the target, without the edges.
    let target_depth_ratio = target.as_ref().and_then(|t| {
        let edge = o.target_tolerance;
        let (lo, hi) = if t.end > t.start + 2 * edge { (t.start + edge, t.end - edge) } else { (t.start, t.end) };
        let depth = c.engine.caller.depth_range(&t.contig, lo, hi)?;
        if depth.is_empty() || genome_depth <= 0.0 { return None }
        Some(depth.iter().map(|&d| d as f64).sum::<f64>() / depth.len() as f64 / genome_depth)
    });

    // The checks.
    let mut checks: Vec<Check> = Vec::new();
    let mut add = |name: &'static str, status: CheckStatus, detail: String| checks.push(Check { name, status, detail });
    let main = sites.first();
    let fmt_site = |s: &Site| format!("{}:{}-{}", s.contig, s.left.map(|v| v.to_string()).unwrap_or_else(|| "?".into()), s.right.map(|v| v.to_string()).unwrap_or_else(|| "?".into()));
    let names = c.constructs.join(", ");
    if sites.is_empty() {
        let status = if enough { CheckStatus::Fail } else { CheckStatus::Inconclusive };
        let detail = if enough { format!("no group of at least {} reads joins the genome to {names}", o.min_junction_reads) } else { format!("no junction was found, and {low}") };
        add("integration_found", status, detail);
    } else {
        add("integration_found", CheckStatus::Pass, format!("{} junction read group(s) join the genome to {names}; the main site is {}", sites.iter().map(|s| (s.left.is_some() as u32) + (s.right.is_some() as u32)).sum::<u32>(), fmt_site(&sites[0])));
    }
    match main {
        None => add("both_sides_joined", CheckStatus::NotApplicable, "no integration site".into()),
        Some(s) if s.left.is_some() && s.right.is_some() => add("both_sides_joined", CheckStatus::Pass, format!("left flank at {}:{} ({} reads), right flank at {}:{} ({} reads)", s.contig, s.left.unwrap_or(0), s.reads_left, s.contig, s.right.unwrap_or(0), s.reads_right)),
        Some(s) => {
            let which = if s.left.is_some() { "left" } else { "right" };
            add("both_sides_joined", if enough { CheckStatus::Fail } else { CheckStatus::Inconclusive }, format!("only the {which} flank is joined to the cassette ({}){}", fmt_site(s), if enough { String::new() } else { format!("; {low}") }));
        }
    }
    match sites.len() {
        0 => add("single_site", CheckStatus::NotApplicable, "no integration site".into()),
        1 => add("single_site", CheckStatus::Pass, format!("one site: {}", fmt_site(&sites[0]))),
        n => add("single_site", CheckStatus::Fail, format!("{n} sites: {}", sites.iter().map(&fmt_site).collect::<Vec<_>>().join("; "))),
    }
    match &target {
        None => {
            add("at_target", CheckStatus::NotApplicable, "no target was given".into());
            add("target_replaced", CheckStatus::NotApplicable, "no target was given".into());
        }
        Some(t) => {
            if let Some(s) = sites.iter().find(|s| at_target(s)) {
                add("at_target", CheckStatus::Pass, format!("the cassette joins the genome at {} (the target is {}:{}-{}, tolerance {} bases)", fmt_site(s), t.contig, t.start, t.end, o.target_tolerance));
            } else if sites.is_empty() {
                add("at_target", if enough { CheckStatus::Fail } else { CheckStatus::Inconclusive }, format!("no integration site was found{}", if enough { String::new() } else { format!(", and {low}") }));
            } else {
                add("at_target", CheckStatus::Fail, format!("the cassette joins the genome at {}, not at the target {}:{}-{}", sites.iter().map(&fmt_site).collect::<Vec<_>>().join("; "), t.contig, t.start, t.end));
            }
            match target_depth_ratio {
                None => add("target_replaced", CheckStatus::Inconclusive, "the depth of the target could not be read".into()),
                Some(_) if !enough => add("target_replaced", CheckStatus::Inconclusive, low.clone()),
                Some(r) if r <= 0.1 => add("target_replaced", CheckStatus::Pass, format!("the inside of the target has {:.0}% of the genome coverage (no wild-type copy)", r * 100.0)),
                Some(r) => add("target_replaced", CheckStatus::Fail, format!("the inside of the target still has {:.0}% of the genome coverage: wild-type sequence is present (a mixed population, an extra copy of the gene, or no replacement)", r * 100.0)),
            }
        }
    }
    if !enough {
        for name in ["cassette_coverage", "copy_number"] {
            add(name, CheckStatus::Inconclusive, low.clone());
        }
        add("cassette_sequence", CheckStatus::Inconclusive, low.clone());
        add("cassette_structure", CheckStatus::Inconclusive, low.clone());
    } else if !any_cassette_reads {
        add("cassette_coverage", CheckStatus::Fail, format!("no reads map to {names}"));
        add("cassette_sequence", CheckStatus::NotApplicable, "no reads map to the cassette".into());
        add("cassette_structure", CheckStatus::NotApplicable, "no reads map to the cassette".into());
        add("copy_number", CheckStatus::Fail, format!("no reads map to {names}: the cassette is not in this sample"));
    } else {
        let gaps: Vec<String> = cassettes.iter().flat_map(|k| k.gaps.iter().map(move |(s, e)| format!("{}:{}-{}", k.name, s, e))).collect();
        if gaps.is_empty() {
            add("cassette_coverage", CheckStatus::Pass, format!("every position of the cassette has at least {gap_depth:.1} reads (genome median {genome_depth:.0}x)"));
        } else {
            add("cassette_coverage", CheckStatus::Fail, format!("positions with fewer than {gap_depth:.1} reads: {}", gaps.join(", ")));
        }
        let variants: Vec<String> = cassettes.iter().flat_map(|k| k.variants.iter().map(move |v| format!("{}:{} {}>{} ({}/{} reads)", k.name, v.pos, v.reference, v.alt, v.alt_count, v.depth))).collect();
        if variants.is_empty() {
            add("cassette_sequence", CheckStatus::Pass, "no difference from the given cassette sequence was called".into());
        } else {
            add("cassette_sequence", CheckStatus::Fail, format!("{} difference(s) from the given sequence: {}", variants.len(), variants.iter().take(10).cloned().collect::<Vec<_>>().join("; ")));
        }
        if internal.is_empty() {
            add("cassette_structure", CheckStatus::Pass, "no read joins one part of the cassette to another".into());
        } else {
            add("cassette_structure", CheckStatus::Fail, format!("reads join parts of the cassette to each other (a deletion, a duplication or a concatemer): {}", internal.iter().map(|j| format!("{}:{} to {}:{} ({} reads)", j.a.contig, j.a.pos, j.b.contig, j.b.pos, j.reads)).collect::<Vec<_>>().join("; ")));
        }
        let ratio = cassettes.first().map(|k| k.copy_ratio).unwrap_or(0.0);
        if (0.6..=1.5).contains(&ratio) {
            add("copy_number", CheckStatus::Pass, format!("cassette depth is {ratio:.2} times the genome depth: one copy per genome"));
        } else if ratio > 1.5 {
            add("copy_number", CheckStatus::Fail, format!("cassette depth is {ratio:.2} times the genome depth: about {:.0} copies per genome", ratio.round().max(2.0)));
        } else {
            add("copy_number", CheckStatus::Fail, format!("cassette depth is only {ratio:.2} times the genome depth: the cassette is in a part of the population only"));
        }
    }
    Report { reads: c.engine.read_count(), genome_depth, junctions, sites, cassettes, checks, target, warnings: c.warnings.clone(), target_depth_ratio }
}

impl Report {
    pub fn check(&self, name: &str) -> Option<&Check> {
        self.checks.iter().find(|c| c.name == name)
    }

    /// Confirmed when no applicable check failed and none was inconclusive.
    pub fn verdict(&self) -> Verdict {
        if self.checks.iter().any(|c| c.status == CheckStatus::Fail) {
            Verdict::NotConfirmed
        } else if self.checks.iter().any(|c| c.status == CheckStatus::Inconclusive) {
            Verdict::Inconclusive
        } else {
            Verdict::Confirmed
        }
    }

    pub fn to_text(&self) -> String {
        let verdict = match self.verdict() {
            Verdict::Confirmed => "CONFIRMED",
            Verdict::NotConfirmed => "NOT CONFIRMED",
            Verdict::Inconclusive => "INCONCLUSIVE",
        };
        let mut out = format!("Allele confirmation: {verdict}\nreads {}, genome median depth {:.1}x\n", self.reads, self.genome_depth);
        if let Some(t) = &self.target {
            out.push_str(&format!("target {}:{}-{}\n", t.contig, t.start, t.end));
        }
        for w in &self.warnings {
            out.push_str(&format!("Warning: {w}\n"));
        }
        out.push_str("\nSites where the cassette joins the genome:\n");
        if self.sites.is_empty() {
            out.push_str("  none\n");
        }
        for s in &self.sites {
            let orient = match s.orientation { Some(Orientation::Forward) => "forward", Some(Orientation::Reverse) => "reverse", None => "unknown" };
            out.push_str(&format!(
                "  {} {}: left flank {} ({} reads, cassette position {}), right flank {} ({} reads, cassette position {}), orientation {orient}, {} genome bases replaced{}\n",
                s.contig, s.construct,
                s.left.map(|v| v.to_string()).unwrap_or_else(|| "not found".into()), s.reads_left, s.construct_at_left.map(|v| v.to_string()).unwrap_or_else(|| "-".into()),
                s.right.map(|v| v.to_string()).unwrap_or_else(|| "not found".into()), s.reads_right, s.construct_at_right.map(|v| v.to_string()).unwrap_or_else(|| "-".into()),
                s.replaced_bases.map(|v| v.to_string()).unwrap_or_else(|| "?".into()),
                if s.homology_left > 0 || s.homology_right > 0 { format!(" (shared bases at the junctions: {} left, {} right; the breakpoints could be that far out)", s.homology_left.max(0), s.homology_right.max(0)) } else { String::new() },
            ));
        }
        out.push_str("\nCassette:\n");
        for k in &self.cassettes {
            out.push_str(&format!("  {} ({} bases): median depth {:.1}, {:.2} times the genome depth, {} coverage gap(s), {} difference(s) from the given sequence\n", k.name, k.length, k.median_depth, k.copy_ratio, k.gaps.len(), k.variants.len()));
        }
        out.push_str("\nChecks:\n");
        for c in &self.checks {
            let tag = match c.status { CheckStatus::Pass => "PASS", CheckStatus::Fail => "FAIL", CheckStatus::Inconclusive => "INCONCLUSIVE", CheckStatus::NotApplicable => "n/a" };
            out.push_str(&format!("  [{tag}] {}: {}\n", c.name, c.detail));
        }
        let constructs: Vec<&str> = self.cassettes.iter().map(|k| k.name.as_str()).collect();
        let side = |s: Side| if s == Side::Left { "left" } else { "right" };
        let line = |j: &Junction| format!("  {}:{} ({}) -- {}:{} ({})  {} reads, overlap {}\n", j.a.contig, j.a.pos, side(j.a.keeps), j.b.contig, j.b.pos, side(j.b.keeps), j.reads, j.overlap);
        let involves = |j: &Junction| constructs.contains(&j.a.contig.as_str()) || constructs.contains(&j.b.contig.as_str());
        let near_target = |j: &Junction| match &self.target {
            Some(t) => [&j.a, &j.b].iter().any(|e| e.contig == t.contig && e.pos + 5_000 >= t.start && e.pos <= t.end + 5_000),
            None => false,
        };
        out.push_str("\nJunctions that involve the cassette (each end is the last reference base next to the breakpoint, and the side of it that the read aligns to):\n");
        let mut any = false;
        for j in self.junctions.iter().filter(|j| involves(j)) {
            out.push_str(&line(j));
            any = true;
        }
        if !any {
            out.push_str("  none\n");
        }
        let near: Vec<&Junction> = self.junctions.iter().filter(|j| !involves(j) && near_target(j)).collect();
        if !near.is_empty() {
            out.push_str("\nJunctions between genome positions within 5 kb of the target:\n");
            near.iter().for_each(|j| out.push_str(&line(j)));
        }
        let others = self.junctions.iter().filter(|j| !involves(j) && !near_target(j)).count();
        if others > 0 {
            out.push_str(&format!("\n{others} other junction(s) between genome positions are in the junctions file only (repeats and mapping artefacts are common there; no check uses them).\n"));
        }
        out
    }

    pub fn to_json(&self) -> String {
        let esc = |s: &str| s.replace('\\', "\\\\").replace('"', "\\\"").replace('\n', "\\n");
        let verdict = match self.verdict() { Verdict::Confirmed => "confirmed", Verdict::NotConfirmed => "not_confirmed", Verdict::Inconclusive => "inconclusive" };
        let status = |s: CheckStatus| match s { CheckStatus::Pass => "pass", CheckStatus::Fail => "fail", CheckStatus::Inconclusive => "inconclusive", CheckStatus::NotApplicable => "not_applicable" };
        let opt = |v: Option<u64>| v.map(|x| x.to_string()).unwrap_or_else(|| "null".to_string());
        let checks: Vec<String> = self.checks.iter().map(|c| format!("{{\"name\":\"{}\",\"status\":\"{}\",\"detail\":\"{}\"}}", c.name, status(c.status), esc(&c.detail))).collect();
        let sites: Vec<String> = self.sites.iter().map(|s| format!(
            "{{\"contig\":\"{}\",\"construct\":\"{}\",\"left\":{},\"right\":{},\"reads_left\":{},\"reads_right\":{},\"orientation\":\"{}\",\"replaced_bases\":{}}}",
            esc(&s.contig), esc(&s.construct), opt(s.left), opt(s.right), s.reads_left, s.reads_right,
            match s.orientation { Some(Orientation::Forward) => "forward", Some(Orientation::Reverse) => "reverse", None => "unknown" },
            s.replaced_bases.map(|v| v.to_string()).unwrap_or_else(|| "null".into()))).collect();
        let warnings: Vec<String> = self.warnings.iter().map(|w| format!("\"{}\"", esc(w))).collect();
        format!("{{\"verdict\":\"{verdict}\",\"reads\":{},\"genome_depth\":{:.2},\"warnings\":[{}],\"sites\":[{}],\"checks\":[{}]}}", self.reads, self.genome_depth, warnings.join(","), sites.join(","), checks.join(","))
    }
}

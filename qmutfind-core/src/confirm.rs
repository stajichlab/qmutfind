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

/// Reads that cross the whole insertion in one piece: a part in the genome on one side, the cassette, and a part in the genome on the other side.
/// A read crosses the whole cassette when it aligns at least this fraction of its bases (long reads lose a few per cent to errors).
const SPANNING_FRACTION: f64 = 0.9;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Spanning {
    /// Reads that join the two flanks of the main site through the cassette.
    pub reads: u32,
    /// Of those, reads that align at least 90% of the cassette bases.
    pub full: u32,
}

/// What the reads say about one construct contig.
#[derive(Debug, Clone, PartialEq)]
pub struct CassetteSummary {
    pub name: String,
    pub length: u64,
    pub median_depth: f64,
    /// Median depth of the cassette divided by the median depth of the genome: about the number of copies per genome.
    pub copy_ratio: f64,
    /// Bases cut off the start and the end of the construct as homology arms (positions on the cassette do not count them).
    pub trimmed_left: u64,
    pub trimmed_right: u64,
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
    /// `Confirmer::from_fasta`: cut homology arms (stretches at the ends of a construct that are also in the genome) off the construct, and take the
    /// target from them when none is given.
    pub strip_arms: bool,
}

impl ConfirmOptions {
    /// The defaults for a kind of reads. Long reads have more errors, so a part of a split read must be longer (100 bases) and junction ends are grouped
    /// within 12 bases; the tolerance between a flank and the target edge is 30 bases.
    pub fn defaults_for(kind: ReadKind) -> Self {
        match kind {
            ReadKind::Short => ConfirmOptions::default(),
            ReadKind::Long => ConfirmOptions { min_part_len: 100, slop: 12, target_tolerance: 30, ..ConfirmOptions::default() },
        }
    }
}

impl Default for ConfirmOptions {
    fn default() -> Self {
        ConfirmOptions { target: None, min_junction_reads: 3, min_part_len: 25, min_part_mapq: 1, slop: 3, target_tolerance: 20, min_depth: 10.0, max_site_span: 200_000, gap_fraction: 0.2, strip_arms: true }
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
    constructs: Vec<String>,
    spanning: Vec<SpanningRead>,
}

/// One read that goes from the genome through a construct and back into the genome.
#[derive(Debug, Clone)]
struct SpanningRead {
    contig: String,
    left: u64,
    right: u64,
    construct: String,
    covered: u64,
}

/// Two parts are neighbours on the read when they overlap by at most this many bases, or are separated by at most this many.
/// A gap on a read, between two aligned parts, of more than `GAP_REALIGN_MIN` and at most `GAP_REALIGN_MAX` bases is aligned again on its own (see `Engine::map_records`).
pub(crate) const GAP_REALIGN_MIN: usize = 200;
pub(crate) const GAP_REALIGN_MAX: usize = 20_000;
const MAX_OVERLAP: i64 = 50;
const MAX_INSERT: i64 = 50;
/// Two parts on one contig, one strand and this close are a gapped alignment, not a junction.
const MIN_SPLIT_DISTANCE: u64 = 50;

impl JunctionCollector {
    pub(crate) fn options(&self) -> JunctionOptions {
        self.options
    }

    pub(crate) fn new(order: Vec<String>, options: JunctionOptions, constructs: Vec<String>) -> Self {
        JunctionCollector { order: order.into_iter().enumerate().map(|(i, n)| (n, i)).collect(), options, observed: Vec::new(), constructs, spanning: Vec::new() }
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
        self.add_spanning(&parts);
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
            // The second part continues the first on the same strand (the read was cut by `Engine::map_records`, or aligned in two pieces): not a junction.
            if a.contig == b.contig && x.reverse == y.reverse && a.pos.abs_diff(b.pos) <= MIN_SPLIT_DISTANCE {
                continue;
            }
            // The same junction from a read of the other strand has its ends in the other order; sort them.
            let key = |e: &End| (self.order.get(&e.contig).copied().unwrap_or(usize::MAX), e.pos, e.keeps);
            if key(&a) > key(&b) {
                std::mem::swap(&mut a, &mut b);
            }
            self.observed.push((a, b, overlap as i32));
        }
    }

    /// Whether two parts that follow each other on the read are neighbours: they overlap by at most `MAX_OVERLAP` bases or are separated by at most `MAX_INSERT`.
    fn neighbours(x: &Alignment, y: &Alignment) -> bool {
        let overlap = x.query_range().1 as i64 - y.query_range().0 as i64;
        overlap <= MAX_OVERLAP && -overlap <= MAX_INSERT
    }

    /// A read with a genome part, then one or more parts of one construct, then a genome part of the same contig, all neighbours on the read.
    fn add_spanning(&mut self, parts: &[&Alignment]) {
        if parts.len() < 3 {
            return;
        }
        let is_construct = |a: &Alignment| self.constructs.iter().any(|c| *c == a.contig);
        for i in 0..parts.len() {
            if is_construct(parts[i]) {
                continue;
            }
            let mut j = i + 1;
            while j < parts.len() && is_construct(parts[j]) && parts[j].contig == parts[i + 1].contig && Self::neighbours(parts[j - 1], parts[j]) {
                j += 1;
            }
            if j > i + 1 && j < parts.len() && !is_construct(parts[j]) && parts[j].contig == parts[i].contig && Self::neighbours(parts[j - 1], parts[j]) {
                let (a, b) = (Self::exit_end(parts[i]), Self::entry_end(parts[j]));
                let middle = &parts[i + 1..j];
                // The cassette bases that the read really aligns to: a stretch that the read skips (a deletion) does not count.
                let covered: u64 = middle.iter().map(|m| m.aligned_ref_bases()).sum();
                self.spanning.push(SpanningRead { contig: a.contig.clone(), left: a.pos.min(b.pos), right: a.pos.max(b.pos), construct: parts[i + 1].contig.clone(), covered });
                return;
            }
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
            // Single linkage: an observation joins the cluster of any earlier one whose two ends are both within `slop` (long reads move a breakpoint by several bases).
            let n = members.len();
            let mut parent: Vec<usize> = (0..n).collect();
            fn find(parent: &mut Vec<usize>, mut i: usize) -> usize {
                while parent[i] != i {
                    parent[i] = parent[parent[i]];
                    i = parent[i];
                }
                i
            }
            for i in 0..n {
                let mut j = i;
                while j > 0 {
                    j -= 1;
                    if members[i].0.pos.abs_diff(members[j].0.pos) > slop {
                        break;
                    }
                    if members[i].1.pos.abs_diff(members[j].1.pos) <= slop {
                        let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                        if ri != rj {
                            parent[ri] = rj;
                        }
                    }
                }
            }
            let mut by_root: HashMap<usize, Vec<&(End, End, i32)>> = HashMap::new();
            for i in 0..n {
                let root = find(&mut parent, i);
                by_root.entry(root).or_default().push(members[i]);
            }
            let mut clusters: Vec<Vec<&(End, End, i32)>> = by_root.into_values().collect();
            clusters.sort_by_key(|c| (c[0].0.pos, c[0].1.pos));
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
/// A shared stretch at the end of a construct is a homology arm when it starts or ends within this many bases of the construct end and is at least `MIN_ARM` long.
const ARM_END_SLACK: usize = 30;
const MIN_ARM: usize = 50;

/// A stretch of a construct (0-based start, end exclusive) that also occurs in the genome, and where: the genome contig and the 1-based interval, and whether
/// the construct runs against the genome strand. `unique` is false when the stretch has several places in the genome.
#[derive(Debug, Clone)]
struct Stretch {
    start: usize,
    end: usize,
    genome: Option<(String, u64, u64, bool)>,
    unique: bool,
}

#[derive(Debug, Clone)]
struct ConstructShare {
    name: String,
    length: usize,
    stretches: Vec<Stretch>,
}

/// For each construct contig: the stretches that are also in the genome (24-base words found in both, either strand), merged, with their genome place.
fn analyze_shared(contigs: &[(String, Vec<u8>)], constructs: &[&str]) -> Vec<ConstructShare> {
    let code = |b: u8| match b.to_ascii_uppercase() { b'A' => Some(0u64), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None };
    let mask = (1u64 << (2 * SHARED_WORD)) - 1;
    // Every word of the constructs, forward and reverse complement: word -> (construct index, start, the word is the reverse complement).
    let mut words: HashMap<u64, Vec<(usize, usize, bool)>> = HashMap::new();
    let mut shares: Vec<(usize, ConstructShare)> = Vec::new();
    for (ci, (name, seq)) in contigs.iter().enumerate() {
        if !constructs.contains(&name.as_str()) { continue }
        shares.push((ci, ConstructShare { name: name.clone(), length: seq.len(), stretches: Vec::new() }));
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
                words.entry(fwd).or_default().push((ci, start, false));
                words.entry(rev).or_default().push((ci, start, true));
            }
        }
    }
    // The genome places of each construct word: (genome contig, 1-based start of the word, the construct word is reversed).
    let mut places: HashMap<(usize, usize), Vec<(usize, u64, bool)>> = HashMap::new();
    let own: Vec<&str> = constructs.to_vec();
    for (gi, (name, seq)) in contigs.iter().enumerate() {
        if own.contains(&name.as_str()) { continue }
        let (mut fwd, mut valid) = (0u64, 0usize);
        for (i, &b) in seq.iter().enumerate() {
            match code(b) {
                Some(c) => { fwd = ((fwd << 2) | c) & mask; valid += 1 }
                None => valid = 0,
            }
            if valid >= SHARED_WORD {
                if let Some(hits) = words.get(&fwd) {
                    for &(ci, start, rev) in hits {
                        let list = places.entry((ci, start)).or_default();
                        if list.len() < 4 { list.push((gi, (i + 2 - SHARED_WORD) as u64, rev)) }
                    }
                }
            }
        }
    }
    let mut out = Vec::new();
    for (ci, mut share) in shares {
        let mut starts: Vec<usize> = places.keys().filter(|(c, _)| *c == ci).map(|&(_, s)| s).collect();
        starts.sort_unstable();
        let mut groups: Vec<(usize, usize, Vec<usize>)> = Vec::new(); // (start, end, word starts)
        for s in starts {
            match groups.last_mut() {
                Some(last) if s <= last.1 => { last.1 = s + SHARED_WORD; last.2.push(s) }
                _ => groups.push((s, s + SHARED_WORD, vec![s])),
            }
        }
        for (start, end, member_starts) in groups {
            // The genome place: the diagonal (genome position minus construct position for the same strand, their sum for the opposite strand) with most words.
            let mut tally: HashMap<(usize, bool, i64), Vec<u64>> = HashMap::new();
            let mut multi = false;
            for s in &member_starts {
                let list = &places[&(ci, *s)];
                if list.len() > 1 { multi = true }
                for &(gi, g, rev) in list {
                    let diag = if rev { g as i64 + *s as i64 } else { g as i64 - *s as i64 };
                    tally.entry((gi, rev, diag)).or_default().push(g);
                }
            }
            let best = tally.iter().max_by_key(|(k, v)| (v.len(), std::cmp::Reverse(**k))).map(|((gi, rev, _), v)| (*gi, *rev, v.clone()));
            let genome = best.map(|(gi, rev, positions)| {
                let lo = *positions.iter().min().unwrap_or(&0);
                let hi = positions.iter().max().copied().unwrap_or(0) + SHARED_WORD as u64 - 1;
                (contigs[gi].0.clone(), lo, hi, rev)
            });
            share.stretches.push(Stretch { start, end, genome, unique: !multi && tally.len() == 1 });
        }
        out.push(share);
    }
    out
}

/// The warnings for constructs that share sequence with the genome: reads from such a stretch map to the genome and to the construct equally well, so the
/// coverage, junctions and copy number of the cassette cannot be trusted. (A construct that was cut by `from_fasta` is analyzed again afterwards, so what
/// is left of such a stretch is still reported.)
fn share_warnings(shares: &[ConstructShare]) -> Vec<String> {
    let mut warnings = Vec::new();
    for share in shares {
        let total: usize = share.stretches.iter().map(|s| s.end - s.start).sum();
        if total >= 50 {
            let list: Vec<String> = share.stretches.iter().filter(|s| s.end - s.start >= 30).take(6).map(|s| format!("{}-{}", s.start + 1, s.end)).collect();
            warnings.push(format!(
                "{total} bases of the construct {} also occur in the genome (positions {}). If these are homology arms, reads from them map to the genome and to the construct equally well, so the coverage, junctions and copy number of the cassette cannot be trusted. Give the cassette sequence without the arms.",
                share.name, list.join(", ")
            ));
        }
    }
    warnings
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
    /// What was done to the construct (arms cut off, target inferred), and a target that came from the arms.
    notes: Vec<String>,
    trimmed: HashMap<String, (u64, u64)>,
    target_inferred: bool,
    kind: ReadKind,
}

impl Confirmer {
    /// From a genome FASTA and a construct FASTA (plain bytes). The construct contigs must not have the name of a genome contig. With `options.strip_arms`,
    /// a stretch of at least 50 bases at the start or the end of a construct that is also in the genome (a homology arm) is cut off; when both ends are
    /// arms, they lie in the same genome contig in the same direction and no target was given, the region between them is the target.
    pub fn from_fasta(genome_fasta: &[u8], construct_fasta: &[u8], kind: ReadKind, mut options: ConfirmOptions) -> Result<Confirmer, String> {
        let mut genome = crate::read_fasta(genome_fasta).map_err(|e| format!("reference FASTA: {e}"))?;
        let mut construct = crate::read_fasta(construct_fasta).map_err(|e| format!("construct FASTA: {e}"))?;
        if construct.is_empty() {
            return Err("the construct FASTA has no sequence: give the sequence you introduced (for example the hph cassette)".to_string());
        }
        for (name, _) in &construct {
            if genome.iter().any(|(g, _)| g == name) {
                return Err(format!("the construct contig {name} has the same name as a contig of the reference; rename it in the construct file"));
            }
        }
        let mut notes = Vec::new();
        let mut trimmed: HashMap<String, (u64, u64)> = HashMap::new();
        let mut inferred = false;
        let mut conflict: Option<String> = None;
        if options.strip_arms {
            let names: Vec<String> = construct.iter().map(|(n, _)| n.clone()).collect();
            let mut with_construct = genome.clone();
            with_construct.extend(construct.iter().cloned());
            let name_refs: Vec<&str> = names.iter().map(String::as_str).collect();
            let shares = analyze_shared(&with_construct, &name_refs);
            drop(with_construct);
            let mut arm_target: Option<Target> = None;
            for share in &shares {
                let left = share.stretches.iter().find(|s| s.start <= ARM_END_SLACK && s.end - s.start >= MIN_ARM);
                let right = share.stretches.iter().rev().find(|s| s.end + ARM_END_SLACK >= share.length && s.end - s.start >= MIN_ARM && left.map_or(true, |l| l.start != s.start));
                let (l, r) = (left.map(|s| s.end).unwrap_or(0), right.map(|s| s.start).unwrap_or(share.length));
                if l >= r {
                    continue; // the whole construct is genome sequence: nothing is cut, and the warning says so
                }
                if left.is_none() && right.is_none() {
                    continue;
                }
                let seq = &mut construct.iter_mut().find(|(n, _)| *n == share.name).expect("the construct is there").1;
                *seq = seq[l..r].to_vec();
                trimmed.insert(share.name.clone(), (l as u64, (share.length - r) as u64));
                let side = |x: Option<&Stretch>| x.map(|s| format!("{} bases", s.end - s.start));
                notes.push(match (side(left), side(right)) {
                    (Some(a), Some(b)) => format!("homology arms were cut off the construct {}: {a} at the start (positions 1-{l}) and {b} at the end (positions {}-{}); positions below are on the cassette without them", share.name, r + 1, share.length),
                    (Some(a), None) => format!("one arm was cut off the construct {}: {a} at the start (positions 1-{l}); no arm was found at the end, so the target is not inferred; positions below are on the cassette without it", share.name),
                    (None, Some(b)) => format!("one arm was cut off the construct {}: {b} at the end (positions {}-{}); no arm was found at the start, so the target is not inferred; positions below are on the cassette without it", share.name, r + 1, share.length),
                    (None, None) => String::new(),
                });
                // The region between the two arms in the genome.
                if let (Some(a), Some(b)) = (left, right) {
                    if let (Some((ca, a1, a2, ra)), Some((cb, b1, b2, rb)), true, true) = (&a.genome, &b.genome, a.unique, b.unique) {
                        if ca == cb && ra == rb {
                            let ((_, first_end), (second_start, _)) = if a1 <= b1 { ((a1, *a2), (*b1, b2)) } else { ((b1, *b2), (*a1, a2)) };
                            if first_end < second_start && second_start - first_end - 1 <= options.max_site_span {
                                arm_target = Some(Target { contig: ca.clone(), start: first_end + 1, end: second_start - 1 });
                            } else if first_end < second_start || (first_end >= second_start && second_start + 1 >= first_end) {
                                notes.push("the two arms are next to each other in the genome: the construct is an insertion, with no region replaced".to_string());
                            }
                        } else {
                            notes.push("the two arms lie in different places or in different directions in the genome, so the target is not inferred".to_string());
                        }
                    } else {
                        notes.push("the arms do not have one clear place in the genome (repeats?), so the target is not inferred".to_string());
                    }
                }
            }
            match (&options.target, arm_target) {
                (None, Some(t)) => {
                    notes.push(format!("the target {}:{}-{} was inferred from the homology arms: it is the genome between them", t.contig, t.start, t.end));
                    options.target = Some(t);
                    inferred = true;
                }
                (Some(given), Some(t)) => {
                    let near = |a: u64, b: u64| a.abs_diff(b) <= options.target_tolerance;
                    if !(given.contig == t.contig && near(given.start, t.start) && near(given.end, t.end)) {
                        conflict = Some(format!(
                            "the homology arms of the construct point to {}:{}-{}, but the target you gave is {}:{}-{}",
                            t.contig, t.start, t.end, given.contig, given.start, given.end
                        ));
                    }
                }
                _ => {}
            }
        }
        let names: Vec<String> = construct.iter().map(|(n, _)| n.clone()).collect();
        genome.extend(construct);
        let mut fasta = Vec::new();
        for (name, seq) in &genome {
            fasta.extend_from_slice(format!(">{name}\n").as_bytes());
            fasta.extend_from_slice(seq);
            fasta.push(b'\n');
        }
        drop(genome);
        let constructs: Vec<&str> = names.iter().map(String::as_str).collect();
        let mut confirmer = Confirmer::new(&fasta, &constructs, kind, options)?;
        confirmer.notes = notes.into_iter().filter(|n| !n.is_empty()).collect();
        confirmer.trimmed = trimmed;
        confirmer.target_inferred = inferred;
        if let Some(c) = conflict {
            confirmer.warnings.insert(0, format!("{c}. The check \"at_target\" uses the target you gave."));
        }
        Ok(confirmer)
    }

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
        engine.track_junctions(JunctionOptions { min_part_len: options.min_part_len, min_part_mapq: options.min_part_mapq }, constructs.iter().map(|s| s.to_string()).collect());
        let contigs = crate::read_fasta(reference_fasta).map_err(|e| format!("reference FASTA: {e}"))?;
        let warnings = share_warnings(&analyze_shared(&contigs, constructs));
        Ok(Confirmer { engine, constructs: constructs.iter().map(|s| s.to_string()).collect(), options, warnings, notes: Vec::new(), trimmed: HashMap::new(), target_inferred: false, kind })
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
    /// Reads that cross the whole insertion (long reads).
    pub spanning: Spanning,
    /// True when the target was not given but taken from the homology arms of the construct.
    pub target_inferred: bool,
    /// What was done to the construct: arms cut off, target inferred.
    pub notes: Vec<String>,
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
    // A one-sided site next to a two-sided one on the same side is the same junction seen with a different breakpoint (long reads): its reads are added to the site.
    let merge = (o.slop * 4).max(40);
    let mut absorbed = vec![false; sites.len()];
    for i in 0..sites.len() {
        if sites[i].left.is_some() && sites[i].right.is_some() {
            continue;
        }
        for k in 0..sites.len() {
            if i == k || absorbed[k] || sites[k].left.is_none() || sites[k].right.is_none() || sites[k].contig != sites[i].contig || sites[k].construct != sites[i].construct {
                continue;
            }
            let near_left = matches!((sites[i].left, sites[k].left), (Some(a), Some(b)) if a.abs_diff(b) <= merge);
            let near_right = matches!((sites[i].right, sites[k].right), (Some(a), Some(b)) if a.abs_diff(b) <= merge);
            if near_left || near_right {
                if near_left { sites[k].reads_left += sites[i].reads_left } else { sites[k].reads_right += sites[i].reads_right }
                absorbed[i] = true;
                break;
            }
        }
    }
    let mut keep = absorbed.into_iter().map(|a| !a);
    sites.retain(|_| keep.next().unwrap_or(true));
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
        let (trimmed_left, trimmed_right) = c.trimmed.get(name).copied().unwrap_or((0, 0));
        cassettes.push(CassetteSummary { name: name.clone(), length, median_depth: median, copy_ratio: if genome_depth > 0.0 { median / genome_depth } else { 0.0 }, trimmed_left, trimmed_right, depth, gaps, variants });
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
                Some(r) if r <= if c.kind == ReadKind::Long { 0.25 } else { 0.1 } => add("target_replaced", CheckStatus::Pass, format!("the inside of the target has {:.0}% of the genome coverage (no wild-type copy)", r * 100.0)),
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
        // Above 1.5 the call depends on the depth noise: short reads at 30x vary by about 20% from place to place, long reads at 15x by more. A second copy
        // measures about 2.0, so a ratio from 1.5 to the failing limit is unclear, never a pass.
        let fail_above = if c.kind == ReadKind::Long { 2.0 } else { 1.75 };
        if (0.6..=1.5).contains(&ratio) {
            add("copy_number", CheckStatus::Pass, format!("cassette depth is {ratio:.2} times the genome depth: one copy per genome"));
        } else if ratio >= fail_above {
            add("copy_number", CheckStatus::Fail, format!("cassette depth is {ratio:.2} times the genome depth: about {:.0} copies per genome", ratio.round().max(2.0)));
        } else if ratio > 1.5 {
            add("copy_number", CheckStatus::Inconclusive, format!("cassette depth is {ratio:.2} times the genome depth: between one and two copies; the depth is too variable at this coverage to tell (a second copy would give about 2)"));
        } else if c.kind == ReadKind::Long && ratio >= 0.35 {
            add("copy_number", CheckStatus::Inconclusive, format!("cassette depth is {ratio:.2} times the genome depth: lower than one copy, but the depth of long reads is too variable at this coverage to tell a mixed population from noise (a half-and-half mixture gives about 0.5)"));
        } else {
            add("copy_number", CheckStatus::Fail, format!("cassette depth is only {ratio:.2} times the genome depth: the cassette is in a part of the population only"));
        }
    }
    // Reads that cross the whole insertion (only long reads can).
    let site_tolerance = o.target_tolerance.max(30);
    let spanning = match sites.first() {
        Some(site) if site.left.is_some() && site.right.is_some() => {
            let (l, r) = (site.left.unwrap_or(0), site.right.unwrap_or(0));
            let cassette_len = cassettes.iter().find(|k| k.name == site.construct).map(|k| k.length).unwrap_or(0);
            let mut sp = Spanning::default();
            for read in &collector.spanning {
                if read.contig == site.contig && read.construct == site.construct && read.left.abs_diff(l) <= site_tolerance && read.right.abs_diff(r) <= site_tolerance {
                    sp.reads += 1;
                    if read.covered as f64 >= SPANNING_FRACTION * cassette_len as f64 {
                        sp.full += 1;
                    }
                }
            }
            sp
        }
        _ => Spanning::default(),
    };
    let two_sided_main = sites.first().map_or(false, |s| s.left.is_some() && s.right.is_some());
    if c.kind == ReadKind::Short {
        add("whole_insertion_in_one_read", CheckStatus::NotApplicable, "short reads (about 150 bases) cannot cross a whole insertion; long reads can show in one piece that the cassette lies between the two flanks".into());
    } else if !two_sided_main {
        add("whole_insertion_in_one_read", CheckStatus::NotApplicable, "there is no site with a junction on both sides".into());
    } else if !enough {
        add("whole_insertion_in_one_read", CheckStatus::Inconclusive, low.clone());
    } else if spanning.full >= 2 {
        add("whole_insertion_in_one_read", CheckStatus::Pass, format!("{} read(s) cross the whole insertion in one piece and align at least 90% of the cassette bases ({} read(s) join the two flanks through it)", spanning.full, spanning.reads));
    } else if spanning.reads >= 1 && spanning.full == 0 {
        add("whole_insertion_in_one_read", CheckStatus::Fail, format!("{} read(s) join the two flanks through the cassette, but none aligns 90% of its bases: part of the cassette is missing", spanning.reads));
    } else if spanning.reads == 0 {
        add("whole_insertion_in_one_read", CheckStatus::Inconclusive, "no read joins the two flanks through the cassette; the reads are probably too short (they must be longer than the cassette plus a part of each flank)".into());
    } else {
        add("whole_insertion_in_one_read", CheckStatus::Inconclusive, format!("only {} read(s) cross the whole insertion; at least 2 are needed", spanning.full));
    }
    Report { reads: c.engine.read_count(), spanning, genome_depth, junctions, sites, cassettes, checks, target, target_inferred: c.target_inferred, notes: c.notes.clone(), warnings: c.warnings.clone(), target_depth_ratio }
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
            out.push_str(&format!("target {}:{}-{}{}\n", t.contig, t.start, t.end, if self.target_inferred { " (inferred from the homology arms of the construct)" } else { "" }));
        }
        for n in &self.notes {
            out.push_str(&format!("Note: {n}\n"));
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
        if self.checks.iter().any(|c| c.name == "whole_insertion_in_one_read" && c.status != CheckStatus::NotApplicable) {
            out.push_str(&format!("\nReads across the whole insertion: {} join the two flanks through the cassette, {} of them align at least 90% of its bases\n", self.spanning.reads, self.spanning.full));
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

    /// All junctions as a table: the two ends, the reads, the overlap, and the kind (`integration`: genome to cassette; `cassette_internal`; `genome`).
    pub fn junctions_tsv(&self) -> String {
        let constructs: Vec<&str> = self.cassettes.iter().map(|k| k.name.as_str()).collect();
        let side = |s: Side| if s == Side::Left { "left" } else { "right" };
        let mut tsv = String::from("contig_a\tpos_a\tside_a\tcontig_b\tpos_b\tside_b\treads\toverlap\tkind\n");
        for j in &self.junctions {
            let kind = match (constructs.contains(&j.a.contig.as_str()), constructs.contains(&j.b.contig.as_str())) {
                (true, true) => "cassette_internal",
                (false, false) => "genome",
                _ => "integration",
            };
            tsv.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{kind}\n", j.a.contig, j.a.pos, side(j.a.keeps), j.b.contig, j.b.pos, side(j.b.keeps), j.reads, j.overlap));
        }
        tsv
    }

    /// The depth at each base of each cassette: contig, position, depth.
    pub fn depth_tsv(&self) -> String {
        let mut out = String::from("contig\tpos\tdepth\n");
        for k in &self.cassettes {
            for (i, d) in k.depth.iter().enumerate() {
                out.push_str(&format!("{}\t{}\t{d}\n", k.name, i + 1));
            }
        }
        out
    }

    pub fn to_json(&self) -> String {
        let esc = |s: &str| s.replace('\\', "\\\\").replace('"', "\\\"").replace('\n', "\\n").replace('\t', "\\t");
        let verdict = match self.verdict() { Verdict::Confirmed => "confirmed", Verdict::NotConfirmed => "not_confirmed", Verdict::Inconclusive => "inconclusive" };
        let status = |s: CheckStatus| match s { CheckStatus::Pass => "pass", CheckStatus::Fail => "fail", CheckStatus::Inconclusive => "inconclusive", CheckStatus::NotApplicable => "not_applicable" };
        let opt = |v: Option<u64>| v.map(|x| x.to_string()).unwrap_or_else(|| "null".to_string());
        let checks: Vec<String> = self.checks.iter().map(|c| format!("{{\"name\":\"{}\",\"status\":\"{}\",\"detail\":\"{}\"}}", c.name, status(c.status), esc(&c.detail))).collect();
        let sites: Vec<String> = self.sites.iter().map(|s| format!(
            "{{\"contig\":\"{}\",\"construct\":\"{}\",\"left\":{},\"right\":{},\"reads_left\":{},\"reads_right\":{},\"orientation\":\"{}\",\"replaced_bases\":{},\"construct_at_left\":{},\"construct_at_right\":{},\"homology_left\":{},\"homology_right\":{}}}",
            esc(&s.contig), esc(&s.construct), opt(s.left), opt(s.right), s.reads_left, s.reads_right,
            match s.orientation { Some(Orientation::Forward) => "forward", Some(Orientation::Reverse) => "reverse", None => "unknown" },
            s.replaced_bases.map(|v| v.to_string()).unwrap_or_else(|| "null".into()), opt(s.construct_at_left), opt(s.construct_at_right), s.homology_left, s.homology_right)).collect();
        let warnings: Vec<String> = self.warnings.iter().map(|w| format!("\"{}\"", esc(w))).collect();
        let cassettes: Vec<String> = self.cassettes.iter().map(|k| {
            let gaps: Vec<String> = k.gaps.iter().map(|(a, b)| format!("[{a},{b}]")).collect();
            let variants: Vec<String> = k.variants.iter().map(|v| format!("{{\"pos\":{},\"ref\":\"{}\",\"alt\":\"{}\",\"alt_count\":{},\"depth\":{}}}", v.pos, esc(&v.reference), esc(&v.alt), v.alt_count, v.depth)).collect();
            format!("{{\"name\":\"{}\",\"length\":{},\"median_depth\":{:.2},\"copy_ratio\":{:.3},\"trimmed_left\":{},\"trimmed_right\":{},\"gaps\":[{}],\"variants\":[{}]}}", esc(&k.name), k.length, k.median_depth, k.copy_ratio, k.trimmed_left, k.trimmed_right, gaps.join(","), variants.join(","))
        }).collect();
        let target = match &self.target {
            Some(t) => format!("{{\"contig\":\"{}\",\"start\":{},\"end\":{}}}", esc(&t.contig), t.start, t.end),
            None => "null".to_string(),
        };
        let notes: Vec<String> = self.notes.iter().map(|w| format!("\"{}\"", esc(w))).collect();
        let ratio = self.target_depth_ratio.map(|r| format!("{r:.3}")).unwrap_or_else(|| "null".into());
        format!(
            "{{\"verdict\":\"{verdict}\",\"reads\":{},\"genome_depth\":{:.2},\"target\":{target},\"spanning\":{{\"reads\":{},\"full\":{}}},\"target_inferred\":{},\"target_depth_ratio\":{ratio},\"notes\":[{}],\"warnings\":[{}],\"sites\":[{}],\"cassettes\":[{}],\"checks\":[{}]}}",
            self.reads, self.genome_depth, self.spanning.reads, self.spanning.full, self.target_inferred, notes.join(","), warnings.join(","), sites.join(","), cassettes.join(","), checks.join(",")
        )
    }
}

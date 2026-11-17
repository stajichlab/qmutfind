//! Breakpoints of an assembly against the reference near one locus: align every assembly contig to the reference (the `rammap`
//! library, a Rust port of minimap2 with its `asm5`, `asm10` and `asm20` presets), keep the alignment blocks near the locus, and
//! report each gap between two blocks of the same contig and strand: the last reference base before it, the first after it, the bases
//! removed from the reference, the bases inserted, and the inserted sequence. The same code runs in the command line tool and in
//! the browser (WebAssembly). The rows are the rows of `breakpoint_summary.tsv`, which `locus_view` draws.

use rammap::{Aligner, Preset, Strand};

use crate::locus::BreakpointRow;
use crate::parse::read_fasta;
use crate::report::Region;

/// What to keep.
#[derive(Debug, Clone, Copy)]
pub struct BreakpointOptions {
    /// Alignment blocks within this many bases of the locus are used.
    pub flank: u64,
    /// Blocks with fewer aligned bases are ignored.
    pub min_aln: u64,
    /// A gap is reported when the reference side or the assembly side is at least this long. An insertion or a deletion of this length inside one
    /// alignment also splits it: the aligner keeps a plain insertion or deletion of a few kb in one alignment, with a long `I` or `D` in its CIGAR.
    pub min_gap: u64,
}

impl Default for BreakpointOptions {
    fn default() -> Self {
        BreakpointOptions { flank: 50_000, min_aln: 1_000, min_gap: 50 }
    }
}

/// The gaps found in one assembly.
#[derive(Debug, Clone)]
pub struct AssemblyBreakpoints {
    pub rows: Vec<BreakpointRow>,
    /// Assembly contigs that have a block of alignment near the locus. Zero means the assembly does not reach the locus.
    pub contigs_near_locus: usize,
}

impl AssemblyBreakpoints {
    /// `{"table":"<breakpoint_summary.tsv text>","fasta":"<inserted sequences>","contigs_near_locus":N}`
    pub fn to_json(&self) -> String {
        format!(
            "{{\"table\":{},\"fasta\":{},\"contigs_near_locus\":{}}}",
            crate::locus::json_string(&breakpoint_table_tsv(&self.rows)),
            crate::locus::json_string(&inserted_sequences_fasta(&self.rows)),
            self.contigs_near_locus
        )
    }
}

/// The index of one reference, to align several assemblies against.
pub struct BreakpointFinder {
    aligner: Aligner,
    contigs: Vec<String>,
    /// The reference sequences (uppercase), to compare the bases of an alignment.
    sequences: Vec<(String, Vec<u8>)>,
}

fn preset_of(name: &str) -> Result<Preset, String> {
    match name {
        "asm5" => Ok(Preset::Asm5),
        "asm10" => Ok(Preset::Asm10),
        "asm20" => Ok(Preset::Asm20),
        other => Err(format!("the preset must be asm5, asm10 or asm20, not {other}")),
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

#[derive(Clone)]
struct Block {
    strand: char,
    qs: usize,
    qe: usize,
    /// 0-based start and exclusive end on the reference; the end is also the 1-based last base.
    ts: usize,
    te: usize,
    mapq: u8,
}

/// Cut one alignment at every stretch where the two sequences do not agree: a long insertion or deletion, or a run of bases that differ. An aligner keeps a plain
/// insertion or deletion of a few kb in one alignment (a long `I` or `D` in the CIGAR), and it aligns a replacement of about the same length straight through
/// as mismatches; only a large difference in length makes it split the alignment. Each alignment column scores -1 for a match and +2 for a mismatch or an
/// inserted or deleted base; a stretch is a run of columns with the highest summed score (the end of the stretch is at its peak), kept when it holds at
/// least `min_gap` bases of the query or of the reference and scores at least `min_gap`.
///
/// `query` and `target` are the aligned stretches (the query already on the strand of the alignment). `ops` are the CIGAR operations (0 match or mismatch,
/// 1 insertion, 2 deletion). Returns the parts between the stretches as `(query start, query end, target start, target end)`: 0-based, end exclusive,
/// query on the forward strand of the contig. A reverse-strand query is read from `query_end` downwards. Without a CIGAR the alignment is one part.
fn split_at_diverged(
    (query_start, query_end): (usize, usize),
    target_start: usize,
    reverse: bool,
    ops: &[(u8, usize)],
    query: &[u8],
    target: &[u8],
    min_gap: usize,
) -> Vec<(usize, usize, usize, usize)> {
    let to_query = |from: usize, to: usize| if reverse { (query_end - to, query_end - from) } else { (query_start + from, query_start + to) };
    // One entry per alignment column: the bad score, and the query and target bases used before it.
    let mut bad: Vec<i8> = Vec::new();
    let mut q_before: Vec<u32> = Vec::new();
    let mut t_before: Vec<u32> = Vec::new();
    let (mut q, mut t) = (0usize, 0usize);
    for &(op, len) in ops {
        for _ in 0..len {
            q_before.push(q as u32);
            t_before.push(t as u32);
            match op {
                0 => {
                    let same = q < query.len() && t < target.len() && query[q] == target[t] && query[q] != b'N';
                    bad.push(if same { -1 } else { 2 });
                    q += 1;
                    t += 1;
                }
                1 => {
                    bad.push(2);
                    q += 1;
                }
                2 => {
                    bad.push(2);
                    t += 1;
                }
                _ => {}
            }
        }
    }
    let columns = bad.len();
    q_before.push(q as u32);
    t_before.push(t as u32);
    // Find the stretches.
    let mut stretches: Vec<(usize, usize)> = Vec::new(); // column ranges, end exclusive
    let (mut cur, mut peak, mut start, mut peak_end) = (0i64, 0i64, 0usize, 0usize);
    let mut close = |start: usize, end: usize, peak: i64, stretches: &mut Vec<(usize, usize)>| {
        let (dq, dt) = ((q_before[end] - q_before[start]) as usize, (t_before[end] - t_before[start]) as usize);
        if peak >= min_gap as i64 && (dq >= min_gap || dt >= min_gap) {
            stretches.push((start, end));
        }
    };
    for (i, &b) in bad.iter().enumerate() {
        cur += b as i64;
        if cur <= 0 {
            if peak > 0 {
                close(start, peak_end, peak, &mut stretches);
            }
            (cur, peak, start, peak_end) = (0, 0, i + 1, i + 1);
        } else if cur > peak {
            (peak, peak_end) = (cur, i + 1);
        } else if peak - cur > 30 {
            // The alignment agrees again: the stretch ended at its peak.
            close(start, peak_end, peak, &mut stretches);
            (cur, peak, start, peak_end) = (0, 0, i + 1, i + 1);
        }
    }
    if peak > 0 {
        close(start, peak_end, peak, &mut stretches);
    }
    // The parts between the stretches.
    let mut parts = Vec::new();
    let mut from = 0usize;
    for &(a, b) in stretches.iter().chain(std::iter::once(&(columns, columns))) {
        if a > from {
            let (qs, qe) = to_query(q_before[from] as usize, q_before[a] as usize);
            let (ts, te) = (target_start + t_before[from] as usize, target_start + t_before[a] as usize);
            if te > ts && qe > qs {
                parts.push((qs, qe, ts, te));
            }
        }
        from = b;
    }
    parts
}

impl BreakpointFinder {
    /// Build the index of `reference_fasta` (plain FASTA bytes). `preset` is `asm5` (assembly and reference less than 5% apart), `asm10` or `asm20`.
    pub fn new(reference_fasta: &[u8], preset: &str) -> Result<BreakpointFinder, String> {
        let reference = read_fasta(reference_fasta).map_err(|e| format!("reference FASTA: {e}"))?;
        if reference.is_empty() {
            return Err("the reference FASTA has no sequences".to_string());
        }
        let contigs = reference.iter().map(|(n, _)| n.clone()).collect();
        let preset = preset_of(preset)?;
        Ok(BreakpointFinder { aligner: Aligner::from_seqs(reference.clone(), preset), contigs, sequences: reference })
    }

    /// Find the gaps of one assembly (plain FASTA bytes) near `locus`. `sample` names the assembly in the rows.
    pub fn assembly(&self, sample: &str, assembly_fasta: &[u8], locus: &Region, options: &BreakpointOptions) -> Result<AssemblyBreakpoints, String> {
        if !self.contigs.iter().any(|c| *c == locus.contig) {
            return Err(format!("the reference has no sequence named {}", locus.contig));
        }
        let contigs = read_fasta(assembly_fasta).map_err(|e| format!("assembly FASTA: {e}"))?;
        if contigs.is_empty() {
            return Err("the assembly FASTA has no sequences".to_string());
        }
        let lo = locus.start.saturating_sub(options.flank);
        let hi = locus.end + options.flank;
        let mut rows = Vec::new();
        let mut near = 0usize;
        for (name, seq) in &contigs {
            let mut blocks: Vec<Block> = Vec::new();
            for m in self.aligner.map_seq(name, seq).mappings.iter() {
                if !m.is_primary || &*m.target_name != locus.contig.as_str() || (m.block_len as u64) < options.min_aln {
                    continue;
                }
                let strand = if matches!(m.strand, Strand::Reverse) { '-' } else { '+' };
                let mapq = m.mapq.clamp(0, 255) as u8;
                let ops: Vec<(u8, usize)> = m.cigar_ops.iter().flatten().map(|c| (c.op as u8, c.len as usize)).collect();
                let target = self.sequences.iter().find(|(n, _)| n == &locus.contig).map(|(_, s)| &s[m.target_start.min(s.len())..m.target_end.min(s.len())]).unwrap_or(&[]);
                let query: Vec<u8> = if strand == '-' { reverse_complement(&seq[m.query_start..m.query_end]) } else { seq[m.query_start..m.query_end].to_vec() };
                for (qs, qe, ts, te) in split_at_diverged((m.query_start, m.query_end), m.target_start, strand == '-', &ops, &query, target, options.min_gap as usize) {
                    if te as u64 >= lo && ts as u64 <= hi {
                        blocks.push(Block { strand, qs, qe, ts, te, mapq });
                    }
                }
            }
            if blocks.is_empty() {
                continue;
            }
            near += 1;
            blocks.sort_by_key(|b| (b.strand, b.ts));
            for pair in blocks.windows(2) {
                let (a, b) = (&pair[0], &pair[1]);
                if a.strand != b.strand {
                    continue;
                }
                let rgap = b.ts as i64 - a.te as i64;
                let qgap = if a.strand == '+' { b.qs as i64 - a.qe as i64 } else { a.qs as i64 - b.qe as i64 };
                if qgap < 0 || rgap < 0 || (qgap == 0 && rgap == 0) || (qgap.max(rgap) as u64) < options.min_gap {
                    continue;
                }
                let (qgap, rgap) = (qgap as u64, rgap as u64);
                let call = if rgap == 0 { "insertion" } else if qgap == 0 { "deletion" } else { "replacement" };
                let (q_lo, q_hi) = if a.strand == '+' { (a.qe, b.qs) } else { (b.qe, a.qs) };
                let ref_left_end = a.te as u64;
                let ref_right_start = b.ts as u64 + 1;
                let (inserted, insert_id) = if qgap > 0 {
                    let piece = &seq[q_lo..q_hi];
                    let piece = if a.strand == '-' { reverse_complement(piece) } else { piece.to_vec() };
                    (
                        Some(String::from_utf8_lossy(&piece).into_owned()),
                        format!("{sample}|{name}:{}-{q_hi}|ref_{ref_left_end}_{ref_right_start}", q_lo + 1),
                    )
                } else {
                    (None, String::new())
                };
                rows.push(BreakpointRow {
                    sample: sample.to_string(),
                    contig: name.clone(),
                    ref_scaffold: locus.contig.clone(),
                    ref_left_end,
                    ref_right_start,
                    ref_gap_bp: rgap,
                    query_gap_bp: qgap,
                    call: call.to_string(),
                    overlaps_locus: ref_left_end < locus.end && ref_right_start > locus.start,
                    strand: a.strand,
                    min_mapq: a.mapq.min(b.mapq),
                    inserted,
                    insert_id,
                });
            }
        }
        rows.sort_by(|a, b| (&a.contig, a.ref_left_end).cmp(&(&b.contig, b.ref_left_end)));
        Ok(AssemblyBreakpoints { rows, contigs_near_locus: near })
    }
}

/// The rows as the text of `breakpoint_summary.tsv`, with a header line.
pub fn breakpoint_table_tsv(rows: &[BreakpointRow]) -> String {
    let mut out = String::from("sample\tcontig\tstrand\tref_scaffold\tref_left_end\tref_right_start\tref_gap_bp\tquery_gap_bp\tcall\toverlaps_locus\tmin_mapq\tinsert_fasta_id\n");
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            r.sample, r.contig, r.strand, r.ref_scaffold, r.ref_left_end, r.ref_right_start, r.ref_gap_bp, r.query_gap_bp, r.call,
            if r.overlaps_locus { "yes" } else { "no" }, r.min_mapq, r.insert_id
        ));
    }
    out
}

/// The inserted sequences of the rows as FASTA (the rows that inserted bases), on the reference strand.
pub fn inserted_sequences_fasta(rows: &[BreakpointRow]) -> String {
    let mut out = String::new();
    for r in rows {
        if let Some(seq) = &r.inserted {
            out.push_str(&format!(">{}\n{seq}\n", r.insert_id));
        }
    }
    out
}

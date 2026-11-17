//! The locus drawing of an insertion or a replacement allele: the reference with its genes on top, the allele below, and the
//! inserted sequence drawn in place of the reference bases that it replaces. The input is one pair of breakpoints, from `qmutfind confirm`
//! (reads) or from a table of assembly breakpoints (`breakpoint_summary.tsv`, one row per gap between two alignment blocks).
//!
//! Both rows use the same number of pixels per base. The allele row is longer or shorter than the reference row by the inserted bases
//! minus the removed bases, so the picture shows what the change does to the length. Genes that overlap the removed bases are listed
//! with what is lost, and their names link to the gene page when a link base is given. Every region in the result is also given as
//! text (`SCAFFOLD:START-END` and `SCAFFOLD:START..END`) that can be copied into a genome browser or into `--target`.

use crate::annotate::{Annotation, Strand, Transcript};
use crate::confirm::{Orientation, Report, Target};
use crate::report::{draw_transcript, escape, group, Region, Scale, LEFT, PLOT, WIDTH};

/// The longest window that is widened to hold whole genes, in bases. A longer gene is drawn clipped.
const MAX_WINDOW: u64 = 200_000;
const GENE_ROW: f64 = 50.0;

/// Where an insertion lies and what it replaces. Positions are 1-based reference positions on `contig`.
#[derive(Debug, Clone, PartialEq)]
pub struct LocusInput {
    pub contig: String,
    /// The last reference base before the insertion.
    pub left: Option<u64>,
    /// The first reference base after the insertion.
    pub right: Option<u64>,
    /// The length of the inserted sequence.
    pub inserted_bases: Option<u64>,
    /// A name for the inserted sequence (the cassette, or the sample of an assembly).
    pub label: String,
    /// The direction of the cassette (`None` when it does not apply, as for an assembly).
    pub orientation: Option<Orientation>,
    /// The region that was meant to be replaced.
    pub target: Option<Target>,
    /// Bases that the reference and the insertion share at each breakpoint (the breakpoint could be that many bases further out).
    pub shared_left: u64,
    pub shared_right: u64,
}

impl LocusInput {
    /// The reference bases that the insertion replaces (`start`, `end`), when both flanks are known and bases lie between them.
    pub fn removed_span(&self) -> Option<(u64, u64)> {
        match (self.left, self.right) {
            (Some(l), Some(r)) if r > l + 1 => Some((l + 1, r - 1)),
            _ => None,
        }
    }

    /// The reference base just before the point where the insertion starts.
    fn point(&self) -> Option<u64> {
        self.left.or_else(|| self.right.map(|r| r.saturating_sub(1)))
    }
}

/// The input of the main insertion site of a `qmutfind confirm` report. `None` when no cassette was found in the genome.
pub fn locus_input_from_report(report: &Report) -> Option<LocusInput> {
    let site = report.sites.first()?;
    let inserted = report.cassettes.iter().find(|c| c.name == site.construct).map(|c| c.length);
    Some(LocusInput {
        contig: site.contig.clone(),
        left: site.left,
        right: site.right,
        inserted_bases: inserted,
        label: site.construct.clone(),
        orientation: site.orientation,
        target: report.target.clone(),
        shared_left: site.homology_left.max(0) as u64,
        shared_right: site.homology_right.max(0) as u64,
    })
}

/// One row of a breakpoint table (`breakpoint_summary.tsv`).
#[derive(Debug, Clone, PartialEq)]
pub struct BreakpointRow {
    pub sample: String,
    /// The assembly contig that carries the insertion.
    pub contig: String,
    pub ref_scaffold: String,
    pub ref_left_end: u64,
    pub ref_right_start: u64,
    pub ref_gap_bp: u64,
    pub query_gap_bp: u64,
    pub call: String,
    /// The table says whether the gap overlaps the locus that the user asked about.
    pub overlaps_locus: bool,
    /// The strand of the alignments: `+` when the assembly contig runs the same way as the reference.
    pub strand: char,
    /// The lower mapping quality of the two alignments that border the gap.
    pub min_mapq: u8,
    /// The inserted sequence on the reference strand. `None` when the row was read from a table, or when nothing was inserted.
    pub inserted: Option<String>,
    /// The name of the inserted sequence in `inserted_sequences.fasta` (empty when nothing was inserted).
    pub insert_id: String,
}

impl BreakpointRow {
    pub fn locus_input(&self) -> LocusInput {
        LocusInput {
            contig: self.ref_scaffold.clone(),
            left: Some(self.ref_left_end),
            right: Some(self.ref_right_start),
            inserted_bases: Some(self.query_gap_bp),
            label: if self.sample.is_empty() { "insertion".to_string() } else { self.sample.clone() },
            orientation: None,
            target: None,
            shared_left: 0,
            shared_right: 0,
        }
    }
}

/// Read a breakpoint table: tab-separated text with a header line. The columns are found by name. Needed: `ref_scaffold`, `ref_left_end`,
/// `ref_right_start`, `ref_gap_bp` and `query_gap_bp`. Used when present: `sample`, `contig`, `call`, `overlaps_locus`.
pub fn parse_breakpoint_table(text: &str) -> Result<Vec<BreakpointRow>, String> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
    let header: Vec<&str> = lines.next().ok_or("the breakpoint table is empty")?.split('\t').map(str::trim).collect();
    let col = |name: &str| header.iter().position(|h| *h == name);
    let need = |name: &str| col(name).ok_or_else(|| format!("the breakpoint table has no column {name}"));
    let (c_scaf, c_left, c_right, c_ref, c_query) = (need("ref_scaffold")?, need("ref_left_end")?, need("ref_right_start")?, need("ref_gap_bp")?, need("query_gap_bp")?);
    let (c_sample, c_contig, c_call, c_overlap) = (col("sample"), col("contig"), col("call"), col("overlaps_locus"));
    let (c_strand, c_mapq, c_id) = (col("strand"), col("min_mapq"), col("insert_fasta_id"));
    let mut rows = Vec::new();
    for (i, line) in lines.enumerate() {
        let f: Vec<&str> = line.split('\t').collect();
        let get = |c: usize| f.get(c).copied().ok_or_else(|| format!("line {} of the breakpoint table has too few columns", i + 2));
        let num = |c: usize, name: &str| get(c)?.trim().parse::<u64>().map_err(|_| format!("line {}: {name} is not a whole number", i + 2));
        let text_of = |c: Option<usize>| c.and_then(|c| f.get(c)).map(|s| s.trim().to_string()).unwrap_or_default();
        rows.push(BreakpointRow {
            sample: text_of(c_sample),
            contig: text_of(c_contig),
            ref_scaffold: get(c_scaf)?.trim().to_string(),
            ref_left_end: num(c_left, "ref_left_end")?,
            ref_right_start: num(c_right, "ref_right_start")?,
            ref_gap_bp: num(c_ref, "ref_gap_bp")?,
            query_gap_bp: num(c_query, "query_gap_bp")?,
            call: text_of(c_call),
            overlaps_locus: text_of(c_overlap) == "yes",
            strand: if text_of(c_strand) == "-" { '-' } else { '+' },
            min_mapq: text_of(c_mapq).parse().unwrap_or(0),
            inserted: None,
            insert_id: text_of(c_id),
        });
    }
    Ok(rows)
}

/// `SCAFFOLD:START-END`, the form that `--target`, the region boxes of the pages and most genome browsers take.
pub fn region_text(contig: &str, start: u64, end: u64) -> String {
    format!("{contig}:{start}-{end}")
}

/// `SCAFFOLD:START..END`, the form of JBrowse and Ensembl.
pub fn region_dots(contig: &str, start: u64, end: u64) -> String {
    format!("{contig}:{start}..{end}")
}

/// A gene that the insertion hits, and what is lost.
#[derive(Debug, Clone, PartialEq)]
pub struct LocusGene {
    pub id: String,
    pub name: Option<String>,
    pub product: Option<String>,
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub gene_bases: u64,
    /// Gene bases inside the removed span (0 when the insertion only interrupts the gene).
    pub removed_bases: u64,
    pub coding_bases: u64,
    pub coding_removed: u64,
    pub effect: String,
}

impl LocusGene {
    pub fn region(&self) -> String {
        region_text(&self.contig, self.start, self.end)
    }
}

/// A named region of the result, to copy.
#[derive(Debug, Clone, PartialEq)]
pub struct NamedRegion {
    pub name: String,
    pub contig: String,
    pub start: u64,
    pub end: u64,
}

impl NamedRegion {
    pub fn text(&self) -> String {
        region_text(&self.contig, self.start, self.end)
    }
    pub fn dots(&self) -> String {
        region_dots(&self.contig, self.start, self.end)
    }
}

/// The drawing, the list of genes that the insertion hits and the regions to copy.
#[derive(Debug, Clone)]
pub struct Locus {
    pub svg: String,
    pub genes: Vec<LocusGene>,
    /// The removed bases, the two flank bases around the insertion (the breakpoints), the target and the drawn window, when they exist.
    pub regions: Vec<NamedRegion>,
}

fn overlap(a: (u64, u64), b: (u64, u64)) -> u64 {
    let (s, e) = (a.0.max(b.0), a.1.min(b.1));
    if e >= s { e - s + 1 } else { 0 }
}

/// What the insertion does to one gene, or `None` when it does not touch it.
fn hit(tx: &Transcript, input: &LocusInput) -> Option<LocusGene> {
    let (gs, ge) = (tx.start(), tx.end());
    let plus = tx.strand != Strand::Minus;
    let coding_bases: u64 = tx.cds.iter().map(|c| c.1 - c.0 + 1).sum();
    let mut gene = LocusGene {
        id: tx.gene.clone(),
        name: tx.gene_name.clone(),
        product: tx.product.clone(),
        contig: tx.contig.clone(),
        start: gs,
        end: ge,
        strand: if plus { '+' } else { '-' },
        gene_bases: ge - gs + 1,
        removed_bases: 0,
        coding_bases,
        coding_removed: 0,
        effect: String::new(),
    };
    if let Some((a, b)) = input.removed_span() {
        gene.removed_bases = overlap((gs, ge), (a, b));
        if gene.removed_bases == 0 {
            return None;
        }
        gene.coding_removed = tx.cds.iter().map(|&c| overlap(c, (a, b))).sum();
        let cuts_left = a <= gs;
        let cuts_right = b >= ge;
        // The 5' end of a gene on the plus strand is its left end.
        let end_name = |left: bool| if left == plus { "5′ end" } else { "3′ end" };
        let mut text = match (cuts_left, cuts_right) {
            (true, true) => "whole gene removed".to_string(),
            _ if coding_bases > 0 && gene.coding_removed == coding_bases => "entire coding sequence removed (UTR sequence outside the span is kept)".to_string(),
            (true, false) => format!("{} removed", end_name(true)),
            (false, true) => format!("{} removed", end_name(false)),
            (false, false) => "internal part removed, gene interrupted".to_string(),
        };
        let whole_cds = coding_bases > 0 && gene.coding_removed == coding_bases;
        if let (false, Some(first), Some(last)) = (whole_cds, tx.cds.first(), tx.cds.last()) {
            // The first coding base in the direction of transcription holds the start codon, the last one the stop codon.
            let (start_base, stop_base) = if plus { (first.0, last.1) } else { (last.1, first.0) };
            if start_base >= a && start_base <= b {
                text.push_str("; start codon removed");
            }
            if stop_base >= a && stop_base <= b {
                text.push_str("; stop codon removed");
            }
        }
        gene.effect = text;
        Some(gene)
    } else {
        let p = input.point()?;
        if gs <= p && p < ge {
            gene.effect = "interrupted by the insertion".to_string();
            Some(gene)
        } else {
            None
        }
    }
}

fn percent_encode(text: &str) -> String {
    let mut out = String::new();
    for b in text.bytes() {
        if b.is_ascii_alphanumeric() || matches!(b, b'-' | b'.' | b'_' | b'~') {
            out.push(b as char);
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}

pub(crate) fn json_string(text: &str) -> String {
    let mut out = String::from("\"");
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

impl Locus {
    /// The genes as a table, with a header line. `region` is `SCAFFOLD:START-END`.
    pub fn genes_tsv(&self) -> String {
        let mut out = String::from("gene\tname\tproduct\tregion\tstrand\tgene_bases\tremoved_bases\tcoding_bases\tcoding_removed\teffect\n");
        for g in &self.genes {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                g.id, g.name.as_deref().unwrap_or(""), g.product.as_deref().unwrap_or(""), g.region(), g.strand, g.gene_bases, g.removed_bases, g.coding_bases, g.coding_removed, g.effect
            ));
        }
        out
    }

    /// The regions as a table: name, `SCAFFOLD:START-END`, `SCAFFOLD:START..END`.
    pub fn regions_tsv(&self) -> String {
        let mut out = String::from("name\tregion\tregion_dots\n");
        for r in &self.regions {
            out.push_str(&format!("{}\t{}\t{}\n", r.name, r.text(), r.dots()));
        }
        out
    }

    /// `{"svg":"...","genes":[...],"regions":[...]}`
    pub fn to_json(&self) -> String {
        let genes: Vec<String> = self
            .genes
            .iter()
            .map(|g| {
                format!(
                    "{{\"id\":{},\"name\":{},\"product\":{},\"region\":{},\"start\":{},\"end\":{},\"strand\":\"{}\",\"gene_bases\":{},\"removed_bases\":{},\"coding_bases\":{},\"coding_removed\":{},\"effect\":{}}}",
                    json_string(&g.id),
                    g.name.as_deref().map(json_string).unwrap_or_else(|| "null".into()),
                    g.product.as_deref().map(json_string).unwrap_or_else(|| "null".into()),
                    json_string(&g.region()),
                    g.start, g.end, g.strand, g.gene_bases, g.removed_bases, g.coding_bases, g.coding_removed,
                    json_string(&g.effect)
                )
            })
            .collect();
        let regions: Vec<String> = self
            .regions
            .iter()
            .map(|r| format!("{{\"name\":{},\"region\":{},\"region_dots\":{}}}", json_string(&r.name), json_string(&r.text()), json_string(&r.dots())))
            .collect();
        format!("{{\"svg\":{},\"genes\":[{}],\"regions\":[{}]}}", json_string(&self.svg), genes.join(","), regions.join(","))
    }
}

/// Draw the locus. `gene_link` is the start of a gene page address (for example `https://fungidb.org/gene/`); the gene ID is added to it.
/// Only genes that the insertion hits get a link.
pub fn locus_view(input: &LocusInput, ann: &Annotation, gene_link: Option<&str>) -> Locus {
    let span = input.removed_span();
    let point = input.point().unwrap_or(1);
    let (core_lo, core_hi) = span.unwrap_or((point, point + 1));

    // The genes that the insertion hits, then the window: the changed bases and the whole of each hit gene, with some space around.
    let mut genes: Vec<LocusGene> = ann
        .transcripts_in(&input.contig, core_lo, core_hi)
        .into_iter()
        .filter_map(|tx| hit(tx, input))
        .collect();
    genes.sort_by_key(|g| (g.start, g.id.clone()));
    let (mut lo, mut hi) = (core_lo, core_hi);
    for g in &genes {
        if g.end.max(hi) - g.start.min(lo) <= MAX_WINDOW {
            lo = lo.min(g.start);
            hi = hi.max(g.end);
        }
    }
    if let Some(t) = input.target.as_ref().filter(|t| t.contig == input.contig) {
        if t.end.max(hi) - t.start.min(lo) <= MAX_WINDOW {
            lo = lo.min(t.start);
            hi = hi.max(t.end);
        }
    }
    let pad = ((hi - lo + 1) / 8).max(300);
    let ws = lo.saturating_sub(pad).max(1);
    let we = hi + pad;
    let region = Region { contig: input.contig.clone(), start: ws, end: we };

    let mut regions = Vec::new();
    let mut add = |name: &str, start: u64, end: u64| regions.push(NamedRegion { name: name.to_string(), contig: input.contig.clone(), start, end });
    if let Some((a, b)) = span {
        add("removed reference bases", a, b);
    }
    if let (Some(l), Some(r)) = (input.left, input.right) {
        add("breakpoints (last base before, first base after)", l, r);
    } else if let Some(p) = input.left.or(input.right) {
        add("breakpoint", p, p);
    }
    if let Some(t) = input.target.as_ref().filter(|t| t.contig == input.contig) {
        add("target", t.start, t.end);
    }
    add("drawn window", ws, we);

    let removed_len = span.map(|(a, b)| b - a + 1).unwrap_or(0);
    let inserted = input.inserted_bases.unwrap_or(0);
    let delta = inserted as i64 - removed_len as i64;
    let ref_len = we - ws + 1;
    let total = ref_len as f64 + delta.max(0) as f64;
    let scale = Scale { start: ws, length: total };
    let px = PLOT / total;

    let transcripts = ann.transcripts_in(&input.contig, ws, we);
    let hit_ids: Vec<&str> = genes.iter().map(|g| g.id.as_str()).collect();
    let mut body = String::new();

    // The axis of the reference row.
    body.push_str(&format!("<line x1=\"{LEFT}\" y1=\"30\" x2=\"{:.1}\" y2=\"30\" style=\"stroke:var(--fg,#222)\"/>", LEFT + ref_len as f64 * px));
    for (pos, anchor) in [(ws, "start"), (ws + (we - ws) / 2, "middle"), (we, "end")] {
        let x = scale.centre(pos);
        body.push_str(&format!(
            "<line x1=\"{x:.1}\" y1=\"26\" x2=\"{x:.1}\" y2=\"34\" style=\"stroke:var(--fg,#222)\"/><text x=\"{x:.1}\" y=\"18\" font-size=\"11\" text-anchor=\"{anchor}\" style=\"fill:var(--fg,#222)\">{}:{}</text>",
            escape(&input.contig),
            group(pos)
        ));
    }
    body.push_str(&format!("<text x=\"{LEFT}\" y=\"48\" font-size=\"11\" style=\"fill:var(--fg,#222)\" font-weight=\"600\">Reference</text>"));

    // The expected target, as a dashed outline behind the genes.
    let genes_top = 56.0;
    let rows = transcripts.len().max(1) as f64;
    let genes_end = genes_top + rows * GENE_ROW;
    if let Some(t) = input.target.as_ref().filter(|t| t.contig == input.contig) {
        let (x1, x2) = scale.span(t.start, t.end);
        body.push_str(&format!(
            "<rect class=\"target\" data-start=\"{}\" data-end=\"{}\" x=\"{x1:.1}\" y=\"{genes_top}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"none\" stroke-dasharray=\"5 3\" style=\"stroke:var(--fg,#222)\" stroke-opacity=\"0.6\"><title>Target {}:{}-{}</title></rect>",
            t.start, t.end, (x2 - x1).max(1.0), genes_end - genes_top, escape(&t.contig), group(t.start), group(t.end)
        ));
    }
    // The removed span, over the genes.
    if let Some((a, b)) = span {
        let (x1, x2) = scale.span(a, b);
        body.push_str(&format!(
            "<rect class=\"removed-span\" data-start=\"{a}\" data-end=\"{b}\" x=\"{x1:.1}\" y=\"{genes_top}\" width=\"{:.1}\" height=\"{:.1}\" style=\"fill:var(--removed,#D55E00)\" fill-opacity=\"0.18\"><title>Removed {} bp: {}:{}-{}</title></rect>",
            (x2 - x1).max(1.0), genes_end - genes_top, group(removed_len), escape(&input.contig), group(a), group(b)
        ));
        body.push_str(&format!(
            "<text x=\"{:.1}\" y=\"48\" font-size=\"11\" text-anchor=\"middle\" style=\"fill:var(--removed,#D55E00)\" font-weight=\"600\">removed {} bp</text>",
            (x1 + x2) / 2.0, group(removed_len)
        ));
    }
    if transcripts.is_empty() {
        body.push_str(&format!("<text x=\"{LEFT}\" y=\"{:.1}\" font-size=\"12\" style=\"fill:var(--fg,#222)\" fill-opacity=\"0.7\">No gene models in this window.</text>", genes_top + 28.0));
    }
    for (i, tx) in transcripts.iter().enumerate() {
        let link = match gene_link {
            Some(base) if hit_ids.contains(&tx.gene.as_str()) => Some(format!("{base}{}", percent_encode(&tx.gene))),
            _ => None,
        };
        draw_transcript(&mut body, tx, &scale, &region, genes_top + i as f64 * GENE_ROW, link.as_deref());
    }

    // The allele row. Left of the insertion it follows the reference; right of it, the reference is shifted by `delta` bases.
    let yb = genes_end + 56.0;
    body.push_str(&format!("<text x=\"{LEFT}\" y=\"{:.1}\" font-size=\"11\" style=\"fill:var(--fg,#222)\" font-weight=\"600\">Allele</text>", yb - 8.0));
    let ins_x0 = LEFT + (point + 1 - ws) as f64 * px;
    let ins_w = if input.inserted_bases.is_some() { (inserted as f64 * px).max(3.0) } else { 40.0 };
    let ins_x1 = ins_x0 + ins_w;
    let flank_style = "style=\"fill:var(--flank,#8a8a8a)\"";
    body.push_str(&format!("<rect class=\"flank-left\" x=\"{LEFT}\" y=\"{yb:.1}\" width=\"{:.1}\" height=\"16\" {flank_style}/>", (ins_x0 - LEFT).max(0.0)));
    let right_known = input.right.is_some() || (span.is_none() && input.left.is_some());
    if right_known {
        let right_first = input.right.unwrap_or(point + 1);
        let width = (we + 1).saturating_sub(right_first) as f64 * px;
        body.push_str(&format!("<rect class=\"flank-right\" x=\"{ins_x1:.1}\" y=\"{yb:.1}\" width=\"{:.1}\" height=\"16\" {flank_style}/>", width.max(0.0)));
    }
    // Lines from the removed span to the inserted block: the block takes the place of the span.
    let (rx1, rx2) = match span {
        Some((a, b)) => scale.span(a, b),
        None => (ins_x0, ins_x0),
    };
    body.push_str(&format!(
        "<polygon class=\"replacement\" points=\"{rx1:.1},{genes_end:.1} {rx2:.1},{genes_end:.1} {ins_x1:.1},{yb:.1} {ins_x0:.1},{yb:.1}\" style=\"fill:var(--cassette,#E69F00)\" fill-opacity=\"0.2\" stroke=\"none\"/>"
    ));
    for (xa, xb) in [(rx1, ins_x0), (rx2, ins_x1)] {
        body.push_str(&format!("<line x1=\"{xa:.1}\" y1=\"{genes_end:.1}\" x2=\"{xb:.1}\" y2=\"{yb:.1}\" style=\"stroke:var(--cassette,#E69F00)\" stroke-width=\"1\"/>"));
    }
    let direction = match input.orientation {
        Some(Orientation::Forward) => "forward",
        Some(Orientation::Reverse) => "reverse",
        None => "",
    };
    let bases_text = if input.inserted_bases.is_some() { format!("{} bp", group(inserted)) } else { "length unknown".to_string() };
    body.push_str(&format!(
        "<rect class=\"insertion\" data-bases=\"{inserted}\" x=\"{ins_x0:.1}\" y=\"{yb:.1}\" width=\"{ins_w:.1}\" height=\"16\" style=\"fill:var(--cassette,#E69F00)\"><title>{} {} inserted{}</title></rect>",
        escape(&input.label), bases_text, if direction.is_empty() { String::new() } else { format!(" ({direction} orientation)") }
    ));
    // Chevrons in the block show the direction of a cassette.
    if let Some(o) = input.orientation {
        let mut cx = ins_x0 + 10.0;
        while cx < ins_x1 - 8.0 {
            let (open, tip) = if o == Orientation::Forward { (cx - 3.0, cx + 3.0) } else { (cx + 3.0, cx - 3.0) };
            body.push_str(&format!("<path class=\"chevron\" d=\"M{open:.1},{:.1} L{tip:.1},{:.1} L{open:.1},{:.1}\" fill=\"none\" stroke=\"#222\" stroke-width=\"1\"/>", yb + 4.0, yb + 8.0, yb + 12.0));
            cx += 22.0;
        }
    }
    // Labels: the two flanks next to the junctions, the inserted sequence below the block.
    let label_y = yb + 32.0;
    let left_text = match input.left {
        Some(l) => format!("left flank ends at {}", group(l)),
        None => "left junction not found".to_string(),
    };
    body.push_str(&format!("<text class=\"junction-left\" x=\"{:.1}\" y=\"{label_y:.1}\" font-size=\"11\" text-anchor=\"end\" style=\"fill:var(--fg,#222)\">{left_text}</text>", ins_x0.max(LEFT + 100.0)));
    let right_text = match input.right {
        Some(r) => format!("right flank starts at {}", group(r)),
        None => "right junction not found".to_string(),
    };
    body.push_str(&format!("<text class=\"junction-right\" x=\"{:.1}\" y=\"{label_y:.1}\" font-size=\"11\" text-anchor=\"start\" style=\"fill:var(--fg,#222)\">{right_text}</text>", ins_x1.min(WIDTH - 180.0)));
    let mut name_line = format!("{} · {}", input.label, bases_text);
    if !direction.is_empty() {
        name_line.push_str(&format!(" · {direction}"));
    }
    body.push_str(&format!(
        "<text class=\"insertion-label\" x=\"{:.1}\" y=\"{:.1}\" font-size=\"12\" text-anchor=\"middle\" style=\"fill:var(--fg,#222)\" font-weight=\"600\">{}</text>",
        ((ins_x0 + ins_x1) / 2.0).clamp(LEFT + 100.0, WIDTH - 100.0), yb + 50.0, escape(&name_line)
    ));
    let mut note_y = yb + 68.0;
    for (shared, side) in [(input.shared_left, "left"), (input.shared_right, "right")] {
        if shared > 0 {
            body.push_str(&format!("<text x=\"{LEFT}\" y=\"{note_y:.1}\" font-size=\"11\" style=\"fill:var(--fg,#222)\" fill-opacity=\"0.8\">The {side} breakpoint is uncertain: the reference and the insertion share {shared} bases there.</text>"));
            note_y += 15.0;
        }
    }
    if span.is_some() || input.inserted_bases.is_some() {
        let net = format!(
            "Removed {} bp, inserted {} bp: the allele is {} bp {} than the reference.",
            group(removed_len), group(inserted), group(delta.unsigned_abs()), if delta < 0 { "shorter" } else { "longer" }
        );
        body.push_str(&format!("<text class=\"net-change\" x=\"{LEFT}\" y=\"{note_y:.1}\" font-size=\"11\" style=\"fill:var(--fg,#222)\">{}</text>", escape(&net)));
        note_y += 15.0;
    }
    let height = note_y + 6.0;
    let svg = format!(
        "<svg class=\"confirm-locus\" xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {WIDTH} {height:.0}\" width=\"100%\" font-family=\"system-ui, sans-serif\" role=\"img\" aria-label=\"Insertion at {}:{}-{}: the reference with its genes, and the allele with the inserted sequence in place of the removed bases\">{body}</svg>",
        escape(&input.contig), ws, we
    );
    Locus { svg, genes, regions }
}

/// Draw the locus with genes read from GFF3 bytes (`None`: no gene models).
pub fn locus_view_gff(input: &LocusInput, gff: Option<&[u8]>, gene_link: Option<&str>) -> Result<Locus, String> {
    let ann = match gff {
        Some(bytes) => Annotation::from_gff(bytes).map_err(|e| format!("GFF: {e}"))?,
        None => Annotation::empty(),
    };
    Ok(locus_view(input, &ann, gene_link))
}

/// The rows of a breakpoint table as JSON, for a menu: `[{"index":1,"sample":..,"contig":..,"ref_scaffold":..,"left":..,"right":..,"ref_gap_bp":..,"query_gap_bp":..,"call":..,"overlaps_locus":true},..]`.
pub fn breakpoint_rows_json(rows: &[BreakpointRow]) -> String {
    let items: Vec<String> = rows
        .iter()
        .enumerate()
        .map(|(i, r)| {
            format!(
                "{{\"index\":{},\"sample\":{},\"contig\":{},\"ref_scaffold\":{},\"left\":{},\"right\":{},\"ref_gap_bp\":{},\"query_gap_bp\":{},\"call\":{},\"overlaps_locus\":{}}}",
                i + 1, json_string(&r.sample), json_string(&r.contig), json_string(&r.ref_scaffold), r.ref_left_end, r.ref_right_start, r.ref_gap_bp, r.query_gap_bp, json_string(&r.call), r.overlaps_locus
            )
        })
        .collect();
    format!("[{}]", items.join(","))
}

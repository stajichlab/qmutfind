//! The locus view: gene models and mutations drawn as SVG inside an HTML page. No dependencies,
//! no JavaScript, one file. Version 1 selects what to draw by a genomic range.

use std::io::BufRead;

use crate::{Annotation, Transcript};

const WIDTH: f64 = 1000.0;
const LEFT: f64 = 60.0;
const PLOT: f64 = 880.0;
/// The width of one character of a label, in pixels at font size 10 (an estimate for the label placement).
const LABEL_CHAR_WIDTH: f64 = 6.0;

/// A genomic range, 1-based and inclusive.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub contig: String,
    pub start: u64,
    pub end: u64,
}

/// Parse `contig:start-end`. Commas in the numbers are allowed.
pub fn parse_region(text: &str) -> Result<Region, String> {
    let bad = || format!("region must look like contig:start-end, not {text}");
    let (contig, range) = text.rsplit_once(':').ok_or_else(bad)?;
    let (start, end) = range.split_once('-').ok_or_else(bad)?;
    let number = |s: &str| s.replace(',', "").parse::<u64>().map_err(|_| bad());
    let (start, end) = (number(start)?, number(end)?);
    if contig.is_empty() || start == 0 || end < start {
        return Err(bad());
    }
    Ok(Region { contig: contig.to_string(), start, end })
}

/// One variant to draw, read from a (usually annotated) VCF.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReportVariant {
    pub contig: String,
    pub pos: u64,
    pub reference: String,
    pub alt: String,
    pub effect: String,
    pub aa: Option<String>,
    pub gene: Option<String>,
    /// The gene's name (for example `hrmA`), when the annotated VCF has one.
    pub gene_name: Option<String>,
    /// The product of the gene (for example `response regulator`).
    pub product: Option<String>,
}

/// Read the variants of a VCF. A record with no `EFFECT` tag is `unannotated`. For a multi-allele
/// record the first allele's effect and protein change are used.
pub fn read_report_variants<R: BufRead>(reader: R) -> Result<Vec<ReportVariant>, String> {
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            return Err(format!("VCF line has {} fields, needs at least 8: {line}", f.len()));
        }
        let info = |key: &str| f[7].split(';').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='));
        out.push(ReportVariant {
            contig: f[0].to_string(),
            pos: f[1].parse().map_err(|_| format!("bad VCF position in: {line}"))?,
            reference: f[3].to_string(),
            alt: f[4].to_string(),
            effect: info("EFFECT").and_then(|e| e.split(',').next()).unwrap_or("unannotated").to_string(),
            aa: info("AA").and_then(|a| a.split(',').next()).filter(|a| *a != ".").map(str::to_string),
            gene: info("GENE").map(str::to_string),
            gene_name: info("GENENAME").map(decode_info),
            product: info("PRODUCT").map(decode_info),
        });
    }
    Ok(out)
}

/// Undo `percent_encode_info`: turn `%XX` back into the character. Anything else is left as it is.
fn decode_info(text: &str) -> String {
    let bytes = text.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 2 < bytes.len() + 0 && bytes[i + 1].is_ascii_hexdigit() && bytes[i + 2].is_ascii_hexdigit() {
            out.push(u8::from_str_radix(&text[i + 1..i + 3], 16).unwrap());
            i += 3;
        } else {
            out.push(bytes[i]);
            i += 1;
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

pub(crate) fn escape(text: &str) -> String {
    text.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;").replace('"', "&quot;")
}

/// A colour-blind-safe palette (Okabe and Ito). Effects without a colour are grey.
/// The short text next to a variant: its protein change (such as `D304G`), else the effect in plain words
/// (`intron`, `frameshift`, `5UTR`). A record with no annotation is `unannotated`.
fn impact_label(v: &ReportVariant) -> String {
    match &v.aa {
        Some(aa) => aa.clone(),
        None => v.effect.trim_end_matches("_derived").replace('_', " "),
    }
}

fn colour(effect: &str) -> &'static str {
    match effect {
        "missense" => "#E69F00",
        "nonsense" | "start_lost" | "stop_lost" => "#D55E00",
        "frameshift" | "large_indel_in_cds" => "#CC79A7",
        "inframe_insertion" | "inframe_deletion" => "#009E73",
        "splice_site" => "#0072B2",
        "5UTR_derived" | "3UTR_derived" | "intron_derived" | "non_coding" => "#56B4E9",
        "synonymous" => "#999999",
        "intergenic" => "#F0E442",
        _ => "#BBBBBB",
    }
}

fn group(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::new();
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

struct Scale {
    start: u64,
    length: f64,
}

impl Scale {
    /// The x position of the centre of base `pos`.
    fn centre(&self, pos: u64) -> f64 {
        LEFT + (pos as f64 - self.start as f64 + 0.5) / self.length * PLOT
    }

    /// The x range covered by bases `from..=to`, clipped to the plot.
    fn span(&self, from: u64, to: u64) -> (f64, f64) {
        let x = |p: f64| (LEFT + (p - self.start as f64) / self.length * PLOT).clamp(LEFT, LEFT + PLOT);
        (x(from as f64), x(to as f64 + 1.0))
    }
}

fn draw_transcript(svg: &mut String, tx: &Transcript, scale: &Scale, region: &Region, y: f64) {
    let strand = if tx.strand == crate::annotate::Strand::Minus { "-" } else { "+" };
    svg.push_str(&format!(
        "<g class=\"transcript\" data-strand=\"{strand}\" data-id=\"{}\">",
        escape(&tx.id)
    ));
    svg.push_str(&format!(
        "<text x=\"{LEFT}\" y=\"{:.1}\" font-size=\"13\" style=\"fill:var(--fg,#222)\">{}</text>",
        y + 6.0,
        escape(&match &tx.gene_name {
            Some(name) => format!("{name} ({})", tx.gene),
            None => tx.gene.clone(),
        })
    ));
    let (x1, x2) = scale.span(tx.exons.first().map(|e| e.0).unwrap_or(0).max(region.start), tx.exons.last().map(|e| e.1).unwrap_or(0).min(region.end));
    svg.push_str(&format!(
        "<line class=\"intron\" x1=\"{x1:.1}\" y1=\"{:.1}\" x2=\"{x2:.1}\" y2=\"{:.1}\" style=\"stroke:var(--fg,#222)\" stroke-width=\"1\"/>",
        y + 22.0,
        y + 22.0
    ));
    let cds_span = tx.cds.first().map(|c| c.0).zip(tx.cds.last().map(|c| c.1));
    for &(s, e) in &tx.exons {
        if e < region.start || s > region.end {
            continue;
        }
        // The exon parts outside the CDS are UTR (derived), or the whole exon of a non-coding transcript.
        let parts: Vec<(u64, u64, &str)> = match cds_span {
            None => vec![(s, e, "noncoding")],
            Some((cs, ce)) => {
                let mut p = Vec::new();
                if s < cs {
                    p.push((s, e.min(cs - 1), "utr"));
                }
                if e > ce {
                    p.push((s.max(ce + 1), e, "utr"));
                }
                p
            }
        };
        for (ps, pe, class) in parts {
            if pe < region.start || ps > region.end {
                continue;
            }
            let (a, b) = scale.span(ps.max(region.start), pe.min(region.end));
            svg.push_str(&format!(
                "<rect class=\"{class}\" x=\"{a:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"8\" style=\"fill:var(--utr,#8aa)\"/>",
                y + 18.0,
                (b - a).max(1.0)
            ));
        }
    }
    for &(s, e) in &tx.cds {
        if e < region.start || s > region.end {
            continue;
        }
        let (a, b) = scale.span(s.max(region.start), e.min(region.end));
        svg.push_str(&format!(
            "<rect class=\"cds\" x=\"{a:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"16\" style=\"fill:var(--cds,#456)\"/>",
            y + 14.0,
            (b - a).max(1.0)
        ));
    }
    svg.push_str("</g>");
}

/// Draw the genes and the variants of one region as SVG. Only variants inside the region are drawn.
pub fn render_locus_svg(ann: &Annotation, variants: &[ReportVariant], region: &Region) -> String {
    let scale = Scale { start: region.start, length: (region.end - region.start + 1) as f64 };
    let transcripts = ann.transcripts_in(&region.contig, region.start, region.end);
    let mut body = String::new();

    // The axis, with the first, middle and last positions.
    body.push_str(&format!("<line x1=\"{LEFT}\" y1=\"30\" x2=\"{:.1}\" y2=\"30\" style=\"stroke:var(--fg,#222)\"/>", LEFT + PLOT));
    for (pos, anchor) in [(region.start, "start"), (region.start + (region.end - region.start) / 2, "middle"), (region.end, "end")] {
        let x = scale.centre(pos);
        body.push_str(&format!(
            "<line x1=\"{x:.1}\" y1=\"26\" x2=\"{x:.1}\" y2=\"34\" style=\"stroke:var(--fg,#222)\"/><text x=\"{x:.1}\" y=\"18\" font-size=\"11\" text-anchor=\"{anchor}\" style=\"fill:var(--fg,#222)\">{}:{}</text>",
            escape(&region.contig),
            group(pos)
        ));
    }

    let mut y = 50.0;
    for tx in &transcripts {
        draw_transcript(&mut body, tx, &scale, region, y);
        y += 50.0;
    }

    // Lollipops. A variant close to the previous one in a level goes up one level.
    let baseline = y + 90.0;
    body.push_str(&format!("<line x1=\"{LEFT}\" y1=\"{baseline}\" x2=\"{:.1}\" y2=\"{baseline}\" style=\"stroke:var(--fg,#222)\" stroke-opacity=\"0.4\"/>", LEFT + PLOT));
    let mut inside: Vec<&ReportVariant> = variants
        .iter()
        .filter(|v| v.contig == region.contig && v.pos >= region.start && v.pos <= region.end)
        .collect();
    inside.sort_by_key(|v| v.pos);
    let mut level_last_x = [f64::NEG_INFINITY; 4];
    let mut used_effects: Vec<String> = Vec::new();
    // First place every circle (x, level), then draw them and the labels. A label needs the places of the others.
    let placed: Vec<(&ReportVariant, f64, usize)> = inside
        .iter()
        .map(|v| {
            let x = scale.centre(v.pos);
            let level = (0..4).find(|&l| x - level_last_x[l] >= 12.0).unwrap_or(3);
            level_last_x[level] = x;
            (*v, x, level)
        })
        .collect();
    let mut label_end = [f64::NEG_INFINITY; 4]; // the right end of the last label drawn on each level
    for (i, (v, x, level)) in placed.iter().enumerate() {
        let (x, level) = (*x, *level);
        let cy = baseline - 14.0 - level as f64 * 16.0;
        let mut tip = format!("{}:{} {}>{} {}", v.contig, v.pos, v.reference, v.alt, v.effect);
        for extra in [&v.aa, &v.gene].into_iter().flatten() {
            tip.push(' ');
            tip.push_str(extra);
        }
        if let Some(name) = &v.gene_name {
            tip.push_str(&format!(" ({name})"));
        }
        body.push_str(&format!(
            "<line x1=\"{x:.1}\" y1=\"{baseline}\" x2=\"{x:.1}\" y2=\"{cy:.1}\" style=\"stroke:var(--fg,#222)\" stroke-opacity=\"0.5\"/>\
<circle class=\"variant\" data-effect=\"{}\" data-pos=\"{}\" cx=\"{x:.2}\" cy=\"{cy:.1}\" r=\"5\" fill=\"{}\" stroke=\"#222\" stroke-width=\"0.5\"><title>{}</title></circle>",
            escape(&v.effect),
            v.pos,
            colour(&v.effect),
            escape(&tip)
        ));
        if !used_effects.contains(&v.effect) {
            used_effects.push(v.effect.clone());
        }

        // The impact label (the protein change, else a short effect name) goes to the right of the circle, or to
        // the left near the right edge. It is left out when it would cover another label or another circle of its
        // level; the tooltip and the table on the page still give it.
        let text = impact_label(v);
        let width = text.chars().count() as f64 * LABEL_CHAR_WIDTH;
        let to_left = x + 8.0 + width > WIDTH - 4.0;
        let (lo, hi) = if to_left { (x - 8.0 - width, x - 8.0) } else { (x + 8.0, x + 8.0 + width) };
        let covers_a_circle = placed.iter().enumerate().any(|(j, (_, x2, l2))| j != i && *l2 == level && *x2 + 5.0 > lo - 2.0 && *x2 - 5.0 < hi + 2.0);
        if lo > label_end[level] + 2.0 && !covers_a_circle {
            label_end[level] = hi;
            body.push_str(&format!(
                "<text class=\"variant-label\" data-pos=\"{}\" x=\"{:.1}\" y=\"{:.1}\" font-size=\"10\" text-anchor=\"{}\" style=\"fill:var(--fg,#222)\">{}</text>",
                v.pos,
                if to_left { x - 8.0 } else { x + 8.0 },
                cy + 3.5,
                if to_left { "end" } else { "start" },
                escape(&text)
            ));
        }
    }

    // The legend lists the effects that appear.
    let legend_y = baseline + 30.0;
    for (i, effect) in used_effects.iter().enumerate() {
        let x = LEFT + (i % 4) as f64 * 220.0;
        let y = legend_y + (i / 4) as f64 * 18.0;
        body.push_str(&format!(
            "<circle class=\"variant-legend\" data-effect=\"{}\" cx=\"{:.1}\" cy=\"{y:.1}\" r=\"5\" fill=\"{}\" stroke=\"#222\" stroke-width=\"0.5\"/><text x=\"{:.1}\" y=\"{:.1}\" font-size=\"12\" style=\"fill:var(--fg,#222)\">{}</text>",
            escape(effect),
            x + 5.0,
            colour(effect),
            x + 16.0,
            y + 4.0,
            escape(effect)
        ));
    }
    let height = legend_y + 18.0 * ((used_effects.len() + 3) / 4) as f64 + 20.0;
    format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {WIDTH} {height:.0}\" width=\"100%\" font-family=\"system-ui, sans-serif\" role=\"img\" aria-label=\"Genes and mutations in {}:{}-{}\">{body}</svg>",
        escape(&region.contig),
        region.start,
        region.end
    )
}

/// Wrap the SVG in a self-contained HTML page with light and dark colours and a summary table.
pub fn render_report_html(title: &str, svg: &str, summary: &[(String, String)]) -> String {
    render_report_html_with_table(title, svg, summary, "")
}

/// The same page with `table_html` (a ready HTML table) after the summary.
pub fn render_report_html_with_table(title: &str, svg: &str, summary: &[(String, String)], table_html: &str) -> String {
    let rows: String = summary
        .iter()
        .map(|(k, v)| format!("<tr><th scope=\"row\">{}</th><td>{}</td></tr>", escape(k), escape(v)))
        .collect();
    format!(
        "<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n<title>{t}</title>\n<style>\n:root {{ --bg: #ffffff; --fg: #1c1c1c; --cds: #3d5a80; --utr: #98a8bc; --line: #d0d0d0; }}\n@media (prefers-color-scheme: dark) {{ :root {{ --bg: #16181d; --fg: #e8e8e8; --cds: #8fb0d8; --utr: #5a6778; --line: #3a3d44; }} }}\nbody {{ background: var(--bg); color: var(--fg); font-family: system-ui, sans-serif; margin: 0 auto; max-width: 1040px; padding: 16px; }}\nh1 {{ font-size: 1.2rem; }}\ntable {{ border-collapse: collapse; margin-top: 12px; }}\nth, td {{ border-bottom: 1px solid var(--line); padding: 4px 12px; text-align: left; }}\n</style></head>\n<body>\n<h1>{t}</h1>\n{svg}\n<table>{rows}</table>\n{table_html}\n</body></html>\n",
        t = escape(title)
    )
}

/// Build the whole page from VCF text, an optional GFF3 (plain bytes) and a region such as
/// `c1:21-110`. The CLI and the browser build both call this, so they give the same page.
pub fn report_page(vcf: &str, gff: Option<&[u8]>, region: &str, title: &str) -> Result<String, String> {
    let region = parse_region(region)?;
    let annotation = match gff {
        Some(bytes) => Annotation::from_gff(bytes).map_err(|e| format!("GFF: {e}"))?,
        None => Annotation::empty(),
    };
    let variants = read_report_variants(vcf.as_bytes())?;
    let in_region = variants
        .iter()
        .filter(|v| v.contig == region.contig && v.pos >= region.start && v.pos <= region.end)
        .count();
    let genes = annotation.transcripts_in(&region.contig, region.start, region.end).len();
    let svg = render_locus_svg(&annotation, &variants, &region);
    let summary = vec![
        ("Region".to_string(), format!("{}:{}-{}", region.contig, group(region.start), group(region.end))),
        ("Variants in the region".to_string(), in_region.to_string()),
        ("Variants in the file".to_string(), variants.len().to_string()),
        ("Genes in the region".to_string(), genes.to_string()),
    ];
    let mut inside: Vec<&ReportVariant> = variants.iter().filter(|v| v.contig == region.contig && v.pos >= region.start && v.pos <= region.end).collect();
    inside.sort_by_key(|v| v.pos);
    let table = if inside.is_empty() {
        String::new()
    } else {
        let rows: String = inside
            .iter()
            .map(|v| {
                format!(
                    "<tr><td>{}:{}</td><td>{} &gt; {}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
                    escape(&v.contig),
                    group(v.pos),
                    escape(&v.reference),
                    escape(&v.alt),
                    escape(&v.effect),
                    escape(v.gene.as_deref().unwrap_or("")),
                    escape(v.gene_name.as_deref().unwrap_or("")),
                    escape(v.aa.as_deref().unwrap_or("")),
                    escape(v.product.as_deref().unwrap_or(""))
                )
            })
            .collect();
        format!("<table class=\"variants\"><thead><tr><th>Position</th><th>Change</th><th>Effect</th><th>Gene</th><th>Name</th><th>Protein change</th><th>Product</th></tr></thead><tbody>{rows}</tbody></table>")
    };
    Ok(render_report_html_with_table(title, &svg, &summary, &table))
}

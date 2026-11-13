//! The whole-genome plot: every call along the contigs, as one SVG. No dependencies, no JavaScript.
//! The top panel counts calls per bin (where the variants are). The bottom panel is a Manhattan-style view with
//! the allele fraction of each call on the y axis. The CLI and the browser build both call this function.

use crate::layout::{bands_and_labels, layout, order_contigs, GAP, LEFT, PLOT};
use crate::report::escape;
use crate::{classify, VariantClass};

const WIDTH: f64 = 1000.0;
const _GAP_USED: f64 = GAP;
const BIN_PX: f64 = 6.0;
const DENSITY_BOTTOM: f64 = 90.0;
const DENSITY_HEIGHT: f64 = 56.0;
const AF_TOP: f64 = 130.0;
const AF_BOTTOM: f64 = 300.0;
const NO_AF_LANE: f64 = 312.0;
const LARGE_INDEL_MIN: usize = 50;

const SNP_COLOUR: &str = "#0072B2";
const SMALL_INDEL_COLOUR: &str = "#E69F00";
const LARGE_INDEL_COLOUR: &str = "#D55E00";

struct Call {
    contig: usize,
    pos: u64,
    reference: String,
    alt: String,
    class: &'static str,
    af: Option<f64>,
}

/// The contigs named in `##contig=<ID=...,length=...>` header lines.
fn header_contigs(vcf: &str) -> Vec<(String, usize)> {
    vcf.lines()
        .filter_map(|l| l.strip_prefix("##contig=<")?.strip_suffix('>'))
        .filter_map(|body| {
            let field = |key: &str| body.split(',').find_map(|kv| kv.strip_prefix(key));
            Some((field("ID=")?.to_string(), field("length=")?.parse().ok()?))
        })
        .collect()
}

fn class_name(reference: &str, alt: &str, tag: Option<&str>) -> &'static str {
    match tag {
        Some("SNP") => "SNP",
        Some("SMALL_INDEL") => "SMALL_INDEL",
        Some("LARGE_INDEL") => "LARGE_INDEL",
        _ if alt.starts_with('<') => "LARGE_INDEL",
        _ => match classify(reference, alt, LARGE_INDEL_MIN) {
            VariantClass::Snp => "SNP",
            VariantClass::SmallIndel => "SMALL_INDEL",
            VariantClass::LargeIndel => "LARGE_INDEL",
        },
    }
}

/// Draw the plot. `contigs` (name and length, in the order to draw) overrides the `##contig` lines of the VCF;
/// pass `None` to use them. A VCF with neither is an error, as is a call on a contig with no known length.
pub fn genome_plot_svg(vcf: &str, contigs: Option<&[(String, usize)]>, order: &str) -> Result<String, String> {
    let contigs: Vec<(String, usize)> = match contigs {
        Some(c) => c.to_vec(),
        None => header_contigs(vcf),
    };
    let contigs = order_contigs(contigs, order)?;
    if contigs.is_empty() {
        return Err("the plot needs contig lengths: the VCF has no ##contig lines, so give the reference genome".to_string());
    }
    let total: usize = contigs.iter().map(|(_, l)| *l).sum();
    if total == 0 {
        return Err("all contigs have length 0".to_string());
    }

    let mut calls: Vec<Call> = Vec::new();
    for line in vcf.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            return Err(format!("VCF line has {} fields, needs at least 8: {line}", f.len()));
        }
        let contig = contigs
            .iter()
            .position(|(n, _)| n == f[0])
            .ok_or_else(|| format!("the VCF has a call on contig {}, which has no length (add it to the contig lines or the reference)", f[0]))?;
        let pos: u64 = f[1].parse().map_err(|_| format!("bad VCF position in: {line}"))?;
        let info = |key: &str| f[7].split(';').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='));
        let alt = f[4].split(',').next().unwrap_or("");
        calls.push(Call {
            contig,
            pos,
            reference: f[3].to_string(),
            alt: alt.to_string(),
            class: class_name(f[3], alt, info("CLASS")),
            af: info("AF").and_then(|a| a.split(',').next()?.parse::<f64>().ok()).map(|a| a.clamp(0.0, 1.0)),
        });
    }

    // Layout: each contig is a band as wide as its share of the genome, with a small gap between bands.
    let l = layout(&contigs);
    let (starts, scale) = (&l.starts, l.scale);
    let band_width = |i: usize| l.widths[i];
    let x_of = |c: &Call| starts[c.contig] + (c.pos.saturating_sub(1)) as f64 * scale;

    // The labels hang below the allele fraction area; the legend and the picture height follow them.
    let label_top = AF_BOTTOM + 10.0;
    let (mut body, label_room) = bands_and_labels(&contigs, &l, 30.0, AF_BOTTOM, label_top);

    // The density panel: calls per bin.
    let mut bins: Vec<Vec<usize>> = contigs.iter().enumerate().map(|(i, _)| vec![0; ((band_width(i) / BIN_PX).round() as usize).max(1)]).collect();
    for c in &calls {
        let len = contigs[c.contig].1.max(1) as f64;
        let n = bins[c.contig].len();
        let idx = (((c.pos.saturating_sub(1)) as f64 / len) * n as f64) as usize;
        bins[c.contig][idx.min(n - 1)] += 1;
    }
    let max_bin = bins.iter().flatten().copied().max().unwrap_or(0);
    for (i, per_contig) in bins.iter().enumerate() {
        let w = band_width(i) / per_contig.len() as f64;
        for (b, &count) in per_contig.iter().enumerate() {
            if count == 0 {
                continue;
            }
            let h = count as f64 / max_bin as f64 * DENSITY_HEIGHT;
            body.push_str(&format!(
                "<rect class=\"bin\" data-count=\"{count}\" x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{h:.2}\" fill=\"#56677f\"><title>{count} call(s)</title></rect>",
                starts[i] + b as f64 * w,
                DENSITY_BOTTOM - h,
                (w - 0.4).max(0.5)
            ));
        }
    }
    body.push_str(&format!("<line x1=\"{LEFT}\" y1=\"{DENSITY_BOTTOM}\" x2=\"{:.1}\" y2=\"{DENSITY_BOTTOM}\" style=\"stroke:var(--fg,#222)\" stroke-opacity=\"0.5\"/>", LEFT + PLOT));
    body.push_str(&format!(
        "<text x=\"{LEFT}\" y=\"22\" font-size=\"12\" style=\"fill:var(--fg,#222)\">Calls per bin (largest bin: {max_bin})</text>"
    ));
    body.push_str(&format!("<text x=\"{LEFT}\" y=\"118\" font-size=\"12\" style=\"fill:var(--fg,#222)\">Allele fraction of each call</text>"));

    // The allele-fraction axis.
    for (af, label) in [(0.0, "0"), (0.5, "0.5"), (1.0, "1")] {
        let y = AF_BOTTOM - af * (AF_BOTTOM - AF_TOP);
        body.push_str(&format!(
            "<line x1=\"{LEFT}\" y1=\"{y:.1}\" x2=\"{:.1}\" y2=\"{y:.1}\" style=\"stroke:var(--fg,#222)\" stroke-opacity=\"{}\"/><text x=\"{:.1}\" y=\"{:.1}\" font-size=\"10\" text-anchor=\"end\" style=\"fill:var(--fg,#222)\">{label}</text>",
            LEFT + PLOT,
            if af == 0.0 { "0.6" } else { "0.15" },
            LEFT - 6.0,
            y + 3.0
        ));
    }

    // The calls.
    let mut lacks_af = false;
    for c in &calls {
        let cx = x_of(c);
        let cy = match c.af {
            Some(af) => AF_BOTTOM - af * (AF_BOTTOM - AF_TOP),
            None => {
                lacks_af = true;
                NO_AF_LANE
            }
        };
        let name = escape(&contigs[c.contig].0);
        let tip = escape(&format!(
            "{}:{} {}>{} {}{}",
            contigs[c.contig].0,
            c.pos,
            c.reference,
            c.alt,
            c.class,
            c.af.map(|a| format!(" AF {a:.2}")).unwrap_or_default()
        ));
        let common = format!("class=\"call\" data-class=\"{}\" data-contig=\"{name}\" data-pos=\"{}\" data-x=\"{cx:.2}\" data-y=\"{cy:.2}\"", c.class, c.pos);
        match c.class {
            "SNP" => body.push_str(&format!("<circle {common} cx=\"{cx:.2}\" cy=\"{cy:.2}\" r=\"3.5\" fill=\"{SNP_COLOUR}\" fill-opacity=\"0.85\" stroke=\"#222\" stroke-width=\"0.4\"><title>{tip}</title></circle>")),
            "SMALL_INDEL" => body.push_str(&format!("<rect {common} x=\"{:.2}\" y=\"{:.2}\" width=\"7\" height=\"7\" fill=\"{SMALL_INDEL_COLOUR}\" fill-opacity=\"0.9\" stroke=\"#222\" stroke-width=\"0.4\"><title>{tip}</title></rect>", cx - 3.5, cy - 3.5)),
            _ => body.push_str(&format!("<polygon {common} points=\"{:.2},{:.2} {:.2},{:.2} {:.2},{:.2}\" fill=\"{LARGE_INDEL_COLOUR}\" fill-opacity=\"0.9\" stroke=\"#222\" stroke-width=\"0.4\"><title>{tip}</title></polygon>", cx, cy - 5.0, cx - 5.0, cy + 4.0, cx + 5.0, cy + 4.0)),
        }
    }
    if calls.is_empty() {
        body.push_str(&format!("<text x=\"{:.1}\" y=\"200\" font-size=\"16\" text-anchor=\"middle\" style=\"fill:var(--fg,#222)\">No variants in this VCF.</text>", LEFT + PLOT / 2.0));
    }
    if lacks_af {
        body.push_str(&format!("<text x=\"{:.1}\" y=\"{:.1}\" font-size=\"10\" text-anchor=\"end\" style=\"fill:var(--fg,#222)\">calls with no AF in the VCF are on this lane</text>", LEFT + PLOT, NO_AF_LANE + 3.0));
    }

    // The legend.
    let ly = label_top + label_room + 22.0;
    body.push_str(&format!("<circle cx=\"{:.1}\" cy=\"{ly}\" r=\"3.5\" fill=\"{SNP_COLOUR}\" stroke=\"#222\" stroke-width=\"0.4\"/><text x=\"{:.1}\" y=\"{:.1}\" font-size=\"11\" style=\"fill:var(--fg,#222)\">SNP</text>", LEFT + 4.0, LEFT + 14.0, ly + 4.0));
    body.push_str(&format!("<rect x=\"{:.1}\" y=\"{:.1}\" width=\"7\" height=\"7\" fill=\"{SMALL_INDEL_COLOUR}\" stroke=\"#222\" stroke-width=\"0.4\"/><text x=\"{:.1}\" y=\"{:.1}\" font-size=\"11\" style=\"fill:var(--fg,#222)\">small indel (1 to 49 bp)</text>", LEFT + 70.0, ly - 3.5, LEFT + 82.0, ly + 4.0));
    body.push_str(&format!("<polygon points=\"{:.1},{:.1} {:.1},{:.1} {:.1},{:.1}\" fill=\"{LARGE_INDEL_COLOUR}\" stroke=\"#222\" stroke-width=\"0.4\"/><text x=\"{:.1}\" y=\"{:.1}\" font-size=\"11\" style=\"fill:var(--fg,#222)\">large indel (50 bp or more)</text>", LEFT + 230.0, ly - 5.0, LEFT + 225.0, ly + 4.0, LEFT + 235.0, ly + 4.0, LEFT + 244.0, ly + 4.0));

    Ok(format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {WIDTH} {:.0}\" width=\"100%\" font-family=\"system-ui, sans-serif\" role=\"img\" aria-label=\"Variants along the genome: {} calls on {} contigs\">{body}</svg>",
        ly + 20.0,
        calls.len(),
        contigs.len()
    ))
}

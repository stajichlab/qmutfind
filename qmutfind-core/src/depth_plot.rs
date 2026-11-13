//! The genome-wide depth plot: mean read depth in bins along the contigs, Manhattan style, as one SVG.
//! The input is BED text with four columns (contig, start, end, depth): the output of `Caller::depth_bed`, or a
//! `mosdepth` regions file, or a `mosdepth` per-base file (rows of any length are averaged by their length).

use crate::layout::{bands_and_labels, layout, order_contigs, LEFT, PLOT};
use crate::report::escape;
use crate::{classify, VariantClass};

const WIDTH: f64 = 1000.0;
const TOP: f64 = 40.0;
const BOTTOM: f64 = 250.0;
const BIN_PX: f64 = 2.0;
const MIN_BIN_BP: f64 = 1000.0;
/// The y axis ends at this multiple of the median depth. Bins above it are drawn at the top and marked.
const Y_MAX_FACTOR: f64 = 2.5;

struct Row {
    contig: usize,
    start: f64,
    end: f64,
    depth: f64,
}

/// Draw the plot. `contigs` (name, length) overrides the lengths found in the BED (the largest end of each contig).
/// `order` is `input`, `length` or `name`. `vcf` adds one tick for each call along the bottom.
pub fn depth_plot_svg(bed: &str, contigs: Option<&[(String, usize)]>, order: &str, vcf: Option<&str>) -> Result<String, String> {
    // Parse the BED.
    let mut parsed: Vec<(String, f64, f64, f64)> = Vec::new();
    for (i, line) in bed.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        let bad = || format!("BED line {}: expected contig, start, end and depth, got: {line}", i + 1);
        if f.len() < 4 {
            return Err(bad());
        }
        let (s, e, d): (f64, f64, f64) = (f[1].parse().map_err(|_| bad())?, f[2].parse().map_err(|_| bad())?, f[3].parse().map_err(|_| bad())?);
        if e <= s || d < 0.0 || !d.is_finite() {
            return Err(bad());
        }
        parsed.push((f[0].to_string(), s, e, d));
    }
    if parsed.is_empty() {
        return Err("the BED has no rows".to_string());
    }
    let contigs: Vec<(String, usize)> = match contigs {
        Some(c) => c.to_vec(),
        None => {
            let mut found: Vec<(String, usize)> = Vec::new();
            for (name, _, e, _) in &parsed {
                match found.iter_mut().find(|(n, _)| n == name) {
                    Some((_, len)) => *len = (*len).max(*e as usize),
                    None => found.push((name.clone(), *e as usize)),
                }
            }
            found
        }
    };
    let contigs = order_contigs(contigs, order)?;
    let total: usize = contigs.iter().map(|(_, l)| *l).sum();
    if total == 0 {
        return Err("all contigs have length 0".to_string());
    }
    let mut rows: Vec<Row> = Vec::with_capacity(parsed.len());
    for (name, s, e, d) in &parsed {
        let contig = contigs.iter().position(|(n, _)| n == name).ok_or_else(|| format!("the BED has a row on contig {name}, which is not in the contig list"))?;
        rows.push(Row { contig, start: *s, end: (*e).min(contigs[contig].1 as f64), depth: *d });
    }

    let l = layout(&contigs);
    // The variant ticks use y 254 to 262, so the labels start below them.
    let label_top = 270.0;
    let (mut body, label_room) = bands_and_labels(&contigs, &l, TOP, BOTTOM, label_top);

    // Bins: about 2 px wide, but never narrower than 1 kb.
    let mut sums: Vec<Vec<f64>> = contigs
        .iter()
        .enumerate()
        .map(|(i, (_, len))| vec![0.0; ((l.widths[i] / BIN_PX).round() as usize).min((*len as f64 / MIN_BIN_BP).ceil() as usize).max(1)])
        .collect();
    for r in &rows {
        let len = contigs[r.contig].1 as f64;
        let n = sums[r.contig].len();
        let bin_bp = len / n as f64;
        let first = ((r.start / bin_bp) as usize).min(n - 1);
        let last = (((r.end - 1e-9).max(r.start) / bin_bp) as usize).min(n - 1);
        for b in first..=last {
            let (lo, hi) = (b as f64 * bin_bp, (b + 1) as f64 * bin_bp);
            let overlap = (r.end.min(hi) - r.start.max(lo)).max(0.0);
            sums[r.contig][b] += r.depth * overlap;
        }
    }
    let means: Vec<Vec<f64>> = sums
        .iter()
        .enumerate()
        .map(|(i, v)| {
            let bin_bp = contigs[i].1 as f64 / v.len() as f64;
            v.iter().map(|s| s / bin_bp).collect()
        })
        .collect();
    let mut all: Vec<f64> = means.iter().flatten().copied().collect();
    all.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if all.len() % 2 == 0 { (all[all.len() / 2 - 1] + all[all.len() / 2]) / 2.0 } else { all[all.len() / 2] };
    let y_max = (median * Y_MAX_FACTOR).max(1.0);
    let y_of = |d: f64| BOTTOM - d.min(y_max) / y_max * (BOTTOM - TOP);

    // The median line and the y axis.
    let my = y_of(median);
    body.push_str(&format!(
        "<line class=\"median\" data-median=\"{median:.2}\" x1=\"{LEFT}\" y1=\"{my:.1}\" x2=\"{:.1}\" y2=\"{my:.1}\" stroke=\"#D55E00\" stroke-dasharray=\"5 3\"/>",
        LEFT + PLOT
    ));
    for (value, label) in [(0.0, "0".to_string()), (median, format!("{median:.0}")), (y_max, format!("{y_max:.0}"))] {
        body.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" font-size=\"10\" text-anchor=\"end\" style=\"fill:var(--fg,#222)\">{label}</text>",
            LEFT - 6.0,
            y_of(value) + 3.0
        ));
    }
    body.push_str(&format!("<line x1=\"{LEFT}\" y1=\"{BOTTOM}\" x2=\"{:.1}\" y2=\"{BOTTOM}\" style=\"stroke:var(--fg,#222)\" stroke-opacity=\"0.6\"/>", LEFT + PLOT));
    body.push_str(&format!(
        "<text x=\"{LEFT}\" y=\"22\" font-size=\"12\" style=\"fill:var(--fg,#222)\">Mean read depth per bin (median {median:.1}x; dashed line)</text>"
    ));

    // One point for each bin.
    for (i, (name, len)) in contigs.iter().enumerate() {
        let n = means[i].len();
        let bin_bp = *len as f64 / n as f64;
        for (b, &d) in means[i].iter().enumerate() {
            let x = l.starts[i] + (b as f64 + 0.5) * bin_bp * l.scale;
            let clipped = d > y_max;
            body.push_str(&format!(
                "<circle class=\"depth\" data-contig=\"{}\" data-depth=\"{d:.2}\"{} data-x=\"{x:.2}\" data-y=\"{:.2}\" cx=\"{x:.2}\" cy=\"{:.2}\" r=\"{}\" fill=\"{}\" fill-opacity=\"0.8\"><title>{}: {:.0}-{:.0} bp, {d:.1}x</title></circle>",
                escape(name),
                if clipped { " data-clipped=\"true\"" } else { "" },
                y_of(d),
                y_of(d),
                if clipped { 3.0 } else { 1.8 },
                if clipped { "#D55E00" } else { "#0072B2" },
                escape(name),
                b as f64 * bin_bp,
                (b + 1) as f64 * bin_bp
            ));
        }
    }

    // Variant ticks along the bottom.
    let mut ticks = 0usize;
    if let Some(vcf) = vcf {
        for line in vcf.lines().filter(|l| !l.starts_with('#') && !l.trim().is_empty()) {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 8 {
                return Err(format!("VCF line has {} fields, needs at least 8: {line}", f.len()));
            }
            let ci = contigs.iter().position(|(n, _)| n == f[0]).ok_or_else(|| format!("the VCF has a call on contig {}, which is not in the contig list", f[0]))?;
            let pos: f64 = f[1].parse().map_err(|_| format!("bad VCF position in: {line}"))?;
            let alt = f[4].split(',').next().unwrap_or("");
            let class = match f[7].split(';').find_map(|kv| kv.strip_prefix("CLASS=")) {
                Some(c) => c.to_string(),
                None if alt.starts_with('<') => "LARGE_INDEL".to_string(),
                None => match classify(f[3], alt, 50) {
                    VariantClass::Snp => "SNP",
                    VariantClass::SmallIndel => "SMALL_INDEL",
                    VariantClass::LargeIndel => "LARGE_INDEL",
                }
                .to_string(),
            };
            let colour = match class.as_str() {
                "SNP" => "#0072B2",
                "SMALL_INDEL" => "#E69F00",
                _ => "#D55E00",
            };
            let x = l.starts[ci] + (pos - 1.0) * l.scale;
            body.push_str(&format!("<line class=\"vtick\" data-class=\"{class}\" x1=\"{x:.2}\" y1=\"254\" x2=\"{x:.2}\" y2=\"262\" stroke=\"{colour}\" stroke-width=\"1\"/>"));
            ticks += 1;
        }
        body.push_str(&format!("<text x=\"{:.1}\" y=\"262\" font-size=\"10\" text-anchor=\"end\" style=\"fill:var(--fg,#222)\">{ticks} variants</text>", LEFT - 6.0));
    }

    Ok(format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {WIDTH} {:.0}\" width=\"100%\" font-family=\"system-ui, sans-serif\" role=\"img\" aria-label=\"Read depth along the genome, {} contigs, median {median:.1}x\">{body}</svg>",
        label_top + label_room + 10.0,
        contigs.len()
    ))
}

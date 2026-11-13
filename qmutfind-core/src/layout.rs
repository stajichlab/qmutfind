//! Shared by the genome-wide plots: how contigs are ordered and laid out side by side along one axis.

use crate::report::escape;

pub(crate) const LEFT: f64 = 70.0;
pub(crate) const PLOT: f64 = 900.0;
pub(crate) const GAP: f64 = 4.0;

/// The x range of each contig and the pixels per base.
pub(crate) struct Layout {
    pub starts: Vec<f64>,
    pub widths: Vec<f64>,
    pub scale: f64,
}

pub(crate) fn layout(contigs: &[(String, usize)]) -> Layout {
    let total: usize = contigs.iter().map(|(_, l)| *l).sum();
    let usable = PLOT - GAP * (contigs.len() as f64 - 1.0);
    let scale = usable / total as f64;
    let mut starts = Vec::with_capacity(contigs.len());
    let mut widths = Vec::with_capacity(contigs.len());
    let mut x = LEFT;
    for (_, len) in contigs {
        starts.push(x);
        let w = *len as f64 * scale;
        widths.push(w);
        x += w + GAP;
    }
    Layout { starts, widths, scale }
}

/// The distance between two rotated labels (the text height plus a gap), the longest label shown, and the width of one bold
/// character at font size 10 (an estimate for the room the labels need), all in pixels or characters.
const LABEL_PITCH: f64 = 12.0;
const MAX_LABEL_CHARS: usize = 24;
const LABEL_CHAR_WIDTH: f64 = 6.6;

/// The alternating shaded bands behind the plot (from `top` to `bottom`) and one label per contig.
///
/// A label is bold and turned 90 degrees (it reads from the bottom up), and right-justified: its end is at (`x`, `label_top`), so
/// the text hangs below `label_top` and never reaches up into the drawing. Two labels are never closer than `LABEL_PITCH`; when
/// two would collide the wider contig keeps its label, and the other shows nothing (its band and label still carry a tooltip with
/// the full name and length). A name over `MAX_LABEL_CHARS` is shortened. Returns the SVG and the height the labels need.
pub(crate) fn bands_and_labels(contigs: &[(String, usize)], l: &Layout, top: f64, bottom: f64, label_top: f64) -> (String, f64) {
    let centres: Vec<f64> = (0..contigs.len()).map(|i| l.starts[i] + l.widths[i] / 2.0).collect();
    // Choose the labels, widest contig first (equal widths: the earlier contig first).
    let mut by_width: Vec<usize> = (0..contigs.len()).collect();
    by_width.sort_by(|&a, &b| l.widths[b].partial_cmp(&l.widths[a]).unwrap_or(std::cmp::Ordering::Equal).then(a.cmp(&b)));
    let mut shown = vec![false; contigs.len()];
    let mut placed: Vec<f64> = Vec::new();
    for i in by_width {
        if placed.iter().all(|c| (c - centres[i]).abs() >= LABEL_PITCH) {
            placed.push(centres[i]);
            shown[i] = true;
        }
    }

    let mut body = String::new();
    let mut longest = 0usize;
    for (i, (name, len)) in contigs.iter().enumerate() {
        body.push_str(&format!(
            "<rect data-band=\"{n}\" x=\"{:.2}\" y=\"{top:.0}\" width=\"{:.2}\" height=\"{:.0}\" style=\"fill:var(--fg,#222)\" fill-opacity=\"{}\"><title>{n}: {len} bp</title></rect>",
            l.starts[i],
            l.widths[i],
            bottom - top,
            if i % 2 == 0 { "0.10" } else { "0.03" },
            n = escape(name)
        ));
        let text: String = if !shown[i] {
            String::new()
        } else if name.chars().count() <= MAX_LABEL_CHARS {
            name.clone()
        } else {
            name.chars().take(MAX_LABEL_CHARS - 1).chain(std::iter::once('…')).collect()
        };
        longest = longest.max(text.chars().count());
        // The text is turned about its anchor; the anchor is moved right by a third of the font size so that the turned text is
        // centred on the contig.
        let x = centres[i] + 3.5;
        body.push_str(&format!(
            "<text data-contig-label=\"{n}\" x=\"{x:.2}\" y=\"{label_top:.0}\" font-size=\"10\" font-weight=\"700\" text-anchor=\"end\" transform=\"rotate(-90 {x:.2} {label_top:.0})\" style=\"fill:var(--fg,#222)\"><title>{n}: {len} bp</title>{}</text>",
            escape(&text),
            n = escape(name)
        ));
    }
    (body, longest as f64 * LABEL_CHAR_WIDTH + 6.0)
}

/// Compare names so that `chr2` comes before `chr10`: digit runs compare as numbers, other text ignores case.
fn natural_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    fn chunks(s: &str) -> Vec<(bool, String)> {
        let mut out: Vec<(bool, String)> = Vec::new();
        for ch in s.chars() {
            let digit = ch.is_ascii_digit();
            match out.last_mut() {
                Some((d, text)) if *d == digit => text.push(ch),
                _ => out.push((digit, ch.to_string())),
            }
        }
        out
    }
    let (ca, cb) = (chunks(a), chunks(b));
    for (x, y) in ca.iter().zip(cb.iter()) {
        let ord = if x.0 && y.0 {
            let (nx, ny) = (x.1.trim_start_matches('0'), y.1.trim_start_matches('0'));
            nx.len().cmp(&ny.len()).then_with(|| nx.cmp(ny))
        } else {
            x.1.to_lowercase().cmp(&y.1.to_lowercase())
        };
        if ord != std::cmp::Ordering::Equal {
            return ord;
        }
    }
    ca.len().cmp(&cb.len()).then_with(|| a.cmp(b))
}

/// Order contigs: `input` keeps the order given, `length` puts the longest first, `name` is natural name order.
pub(crate) fn order_contigs(mut contigs: Vec<(String, usize)>, how: &str) -> Result<Vec<(String, usize)>, String> {
    match how {
        "input" => {}
        "length" => contigs.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| natural_cmp(&a.0, &b.0))),
        "name" => contigs.sort_by(|a, b| natural_cmp(&a.0, &b.0)),
        other => return Err(format!("the contig order must be input, length or name, not {other}")),
    }
    Ok(contigs)
}

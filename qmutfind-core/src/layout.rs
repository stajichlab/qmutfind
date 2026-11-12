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

/// The alternating shaded bands behind the plot (from `top` to `bottom`) and one label per contig at `label_y`.
/// A label that does not fit is shortened; the full name and length stay in a tooltip.
pub(crate) fn bands_and_labels(contigs: &[(String, usize)], l: &Layout, top: f64, bottom: f64, label_y: f64) -> String {
    let mut body = String::new();
    for (i, (name, len)) in contigs.iter().enumerate() {
        body.push_str(&format!(
            "<rect data-band=\"{}\" x=\"{:.2}\" y=\"{top:.0}\" width=\"{:.2}\" height=\"{:.0}\" style=\"fill:var(--fg,#222)\" fill-opacity=\"{}\"/>",
            escape(name),
            l.starts[i],
            l.widths[i],
            bottom - top,
            if i % 2 == 0 { "0.10" } else { "0.03" }
        ));
        let fit = (l.widths[i] / 6.5) as usize;
        let shown: String = if name.chars().count() <= fit {
            name.clone()
        } else if fit >= 3 {
            name.chars().take(fit - 1).chain(std::iter::once('…')).collect()
        } else {
            String::new()
        };
        body.push_str(&format!(
            "<text data-contig-label=\"{n}\" x=\"{:.2}\" y=\"{label_y:.0}\" font-size=\"10\" text-anchor=\"middle\" style=\"fill:var(--fg,#222)\"><title>{n}: {len} bp</title>{}</text>",
            l.starts[i] + l.widths[i] / 2.0,
            escape(&shown),
            n = escape(name)
        ));
    }
    body
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

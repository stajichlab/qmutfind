//! Helpers shared by the integration tests. The reference is synthetic and deterministic.
#![allow(dead_code)]

/// Small deterministic generator, so the reference is the same on every machine.
pub struct Lcg(pub u64);

impl Lcg {
    pub fn base(&mut self) -> u8 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        b"ACGT"[((self.0 >> 33) % 4) as usize]
    }
}

pub fn reference() -> Vec<(String, Vec<u8>)> {
    let mut rng = Lcg(42);
    let chr_a: Vec<u8> = (0..120_000).map(|_| rng.base()).collect();
    let chr_b: Vec<u8> = (0..80_000).map(|_| rng.base()).collect();
    vec![("chrA".to_string(), chr_a), ("chrB".to_string(), chr_b)]
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
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

/// A different base from `b`.
pub fn other_base(b: u8) -> u8 {
    if b == b'G' { b'T' } else { b'G' }
}

// ---- The contig labels of the genome-wide plots --------------------------------------------------------------

pub struct ContigLabel {
    pub name: String,
    pub x: f64,
    pub y: f64,
    pub text: String,
    pub anchor: String,
    pub weight: String,
    pub transform: String,
}

fn attribute(header: &str, name: &str) -> String {
    let key = format!(" {name}=\"");
    let at = header.find(&key).unwrap_or_else(|| panic!("no {name} in {header}")) + key.len();
    header[at..at + header[at..].find('"').unwrap()].to_string()
}

/// Every `<text data-contig-label=...>` of a plot: where it is, how it is drawn, and the text that is shown (without the tooltip).
pub fn contig_labels(svg: &str) -> Vec<ContigLabel> {
    svg.match_indices("<text data-contig-label=\"")
        .map(|(i, _)| {
            let header = &svg[i..i + svg[i..].find('>').unwrap()];
            let body = &svg[i + header.len() + 1..i + svg[i..].find("</text>").unwrap()];
            let text = match body.find("</title>") {
                Some(t) => &body[t + "</title>".len()..],
                None => body,
            };
            ContigLabel {
                name: attribute(header, "data-contig-label"),
                x: attribute(header, "x").parse().unwrap(),
                y: attribute(header, "y").parse().unwrap(),
                text: text.to_string(),
                anchor: attribute(header, "text-anchor"),
                weight: attribute(header, "font-weight"),
                transform: attribute(header, "transform"),
            }
        })
        .collect()
}

pub fn viewbox_height(svg: &str) -> f64 {
    let at = svg.find("viewBox=\"").unwrap() + "viewBox=\"".len();
    svg[at..at + svg[at..].find('"').unwrap()].split_whitespace().nth(3).unwrap().parse().unwrap()
}

/// 55 contigs like a scaffold assembly: a few long ones and many tiny ones.
pub fn scaffold_contigs() -> Vec<(String, usize)> {
    (0..55).map(|i| (format!("DS4995{:02}.1", i), if i < 8 { 3_000_000 } else { 2_000 + i * 700 })).collect()
}

/// A long filler, then a very narrow contig, then a wider one whose label centre is 11 px from the narrow one's.
pub fn crowded_contigs() -> Vec<(String, usize)> {
    vec![("filler".to_string(), 1_000_000), ("narrow".to_string(), 2_300), ("wider".to_string(), 13_600)]
}

/// The rules for the labels of one plot: bold, rotated 90 degrees, right-justified so that they hang below `min_y`, no two
/// labels closer than the text height, and a plot tall enough for the longest label.
pub fn assert_label_rules(svg: &str, min_y: f64) {
    let labels = contig_labels(svg);
    assert!(!labels.is_empty());
    let mut height_needed: f64 = 0.0;
    for l in &labels {
        assert_eq!(l.weight, "700", "{} is not bold", l.name);
        assert_eq!(l.anchor, "end", "{} is not right-justified", l.name);
        let nums: Vec<f64> = l.transform.strip_prefix("rotate(").and_then(|t| t.strip_suffix(')')).unwrap_or("").split_whitespace().map(|n| n.parse().unwrap()).collect();
        assert!(nums.len() == 3 && nums[0] == -90.0 && (nums[1] - l.x).abs() < 0.01 && (nums[2] - l.y).abs() < 0.01, "{} is not turned 90 degrees about its own anchor: {}", l.name, l.transform);
        assert!(l.y >= min_y, "{} starts at y {} above the limit {min_y}", l.name, l.y);
        if !l.text.is_empty() {
            height_needed = height_needed.max(l.y + l.text.chars().count() as f64 * 6.6);
        }
    }
    let shown: Vec<&ContigLabel> = labels.iter().filter(|l| !l.text.is_empty()).collect();
    for (i, a) in shown.iter().enumerate() {
        for b in &shown[i + 1..] {
            assert!((a.x - b.x).abs() >= 11.5, "{} and {} are {} px apart", a.name, b.name, (a.x - b.x).abs());
        }
    }
    assert!(viewbox_height(svg) >= height_needed, "the plot is {} high but the longest label needs {height_needed}", viewbox_height(svg));
}

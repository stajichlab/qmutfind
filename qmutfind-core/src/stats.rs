//! Summary statistics of a variant set, read from a VCF. These need only the VCF (and the contig
//! lengths in its header). Statistics that need the read depth, the reference sequence or per-strand
//! counts are not here.

use std::collections::BTreeMap;
use std::io::BufRead;

use crate::classify;

#[derive(Debug, Clone, Default)]
pub struct ContigStats {
    pub length: u64,
    pub count: u64,
    /// Variants per kb of contig. 0 when the contig length is not in the header.
    pub per_kb: f64,
}

#[derive(Debug, Clone, Default)]
pub struct Stats {
    pub total: u64,
    pub by_class: BTreeMap<String, u64>,
    pub by_effect: BTreeMap<String, u64>,
    pub transitions: u64,
    pub transversions: u64,
    pub ts_tv: Option<f64>,
    /// The 12 substitution types, for example `A>G`.
    pub spectrum: BTreeMap<String, u64>,
    /// The 6 classes with the pyrimidine of each pair as the reference, for example `T>C`.
    pub spectrum_collapsed: BTreeMap<String, u64>,
    /// Signed length (insertion positive, deletion negative) to count.
    pub indel_lengths: BTreeMap<i64, u64>,
    /// Share of indels whose length is not a multiple of three.
    pub frameshift_length_fraction: Option<f64>,
    pub per_contig: BTreeMap<String, ContigStats>,
    pub overall_per_kb: Option<f64>,
    /// Ten bins of alternate fraction, 0.0 to 0.1 up to 0.9 to 1.0 (1.0 is in the last bin).
    pub fraction_histogram: Vec<u64>,
    /// Records with a fraction from 0.2 to 0.8. A haploid sample should have almost none.
    pub mid_fraction_calls: u64,
    pub mid_fraction_share: Option<f64>,
    /// Smallest distance between neighbouring variants on one contig.
    pub smallest_gap: Option<u64>,
    /// Variants with a neighbour within `cluster_window` bases.
    pub clustered: u64,
}

fn is_transition(r: u8, a: u8) -> bool {
    matches!((r, a), (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C'))
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        o => o,
    }
}

/// Compute the statistics. `cluster_window` is the distance in bases for the clustering count.
pub fn compute_stats<R: BufRead>(reader: R, cluster_window: u64) -> Result<Stats, String> {
    let mut s = Stats { fraction_histogram: vec![0; 10], ..Stats::default() };
    let mut positions: BTreeMap<String, Vec<u64>> = BTreeMap::new();
    let mut with_fraction = 0u64;
    let mut indels = 0u64;
    let mut frameshift_lengths = 0u64;

    for line in reader.lines() {
        let line = line.map_err(|e| e.to_string())?;
        if let Some(rest) = line.strip_prefix("##contig=<") {
            let attr = |key: &str| rest.trim_end_matches('>').split(',').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('=').map(str::to_string));
            if let (Some(id), Some(len)) = (attr("ID"), attr("length").and_then(|l| l.parse().ok())) {
                s.per_contig.entry(id).or_default().length = len;
            }
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            return Err(format!("VCF line has {} fields, needs at least 8: {line}", f.len()));
        }
        let pos: u64 = f[1].parse().map_err(|_| format!("bad VCF position in: {line}"))?;
        let info = |key: &str| f[7].split(';').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='));
        let (reference, alt) = (f[3], f[4].split(',').next().unwrap_or(""));

        s.total += 1;
        let class = match info("CLASS") {
            Some(c) => c.to_string(),
            None => match classify(reference, alt, 50) {
                crate::VariantClass::Snp => "SNP",
                crate::VariantClass::SmallIndel => "SMALL_INDEL",
                crate::VariantClass::LargeIndel => "LARGE_INDEL",
            }
            .to_string(),
        };
        *s.by_class.entry(class).or_default() += 1;
        *s.by_effect.entry(info("EFFECT").unwrap_or("unannotated").to_string()).or_default() += 1;
        s.per_contig.entry(f[0].to_string()).or_default().count += 1;
        positions.entry(f[0].to_string()).or_default().push(pos);

        if reference.len() == 1 && alt.len() == 1 {
            let (r, a) = (reference.as_bytes()[0].to_ascii_uppercase(), alt.as_bytes()[0].to_ascii_uppercase());
            if matches!(r, b'A' | b'C' | b'G' | b'T') && matches!(a, b'A' | b'C' | b'G' | b'T') && r != a {
                if is_transition(r, a) {
                    s.transitions += 1;
                } else {
                    s.transversions += 1;
                }
                *s.spectrum.entry(format!("{}>{}", r as char, a as char)).or_default() += 1;
                let (cr, ca) = if matches!(r, b'A' | b'G') { (complement(r), complement(a)) } else { (r, a) };
                *s.spectrum_collapsed.entry(format!("{}>{}", cr as char, ca as char)).or_default() += 1;
            }
        } else {
            let length = if alt.starts_with('<') {
                info("SVLEN").and_then(|v| v.split(',').next()?.parse::<i64>().ok())
            } else {
                Some(alt.len() as i64 - reference.len() as i64)
            };
            if let Some(length) = length.filter(|l| *l != 0) {
                *s.indel_lengths.entry(length).or_default() += 1;
                indels += 1;
                if length % 3 != 0 {
                    frameshift_lengths += 1;
                }
            }
        }

        if let Some(af) = info("AF").and_then(|v| v.split(',').next()?.parse::<f64>().ok()) {
            with_fraction += 1;
            // Integer arithmetic on the fraction rounded to four decimals, so 0.7 cannot land in the wrong bin.
            let milli = (af * 10_000.0).round() as i64;
            let bin = ((milli / 1000).clamp(0, 9)) as usize;
            s.fraction_histogram[bin] += 1;
            if (2000..=8000).contains(&milli) {
                s.mid_fraction_calls += 1;
            }
        }
    }

    s.ts_tv = if s.transversions > 0 { Some(s.transitions as f64 / s.transversions as f64) } else { None };
    s.frameshift_length_fraction = if indels > 0 { Some(frameshift_lengths as f64 / indels as f64) } else { None };
    s.mid_fraction_share = if with_fraction > 0 { Some(s.mid_fraction_calls as f64 / with_fraction as f64) } else { None };
    let total_length: u64 = s.per_contig.values().map(|c| c.length).sum();
    s.overall_per_kb = if total_length > 0 { Some(s.total as f64 / (total_length as f64 / 1000.0)) } else { None };
    for c in s.per_contig.values_mut() {
        if c.length > 0 {
            c.per_kb = c.count as f64 / (c.length as f64 / 1000.0);
        }
    }
    for list in positions.values_mut() {
        list.sort_unstable();
        let gaps: Vec<u64> = list.windows(2).map(|w| w[1] - w[0]).collect();
        s.smallest_gap = match (s.smallest_gap, gaps.iter().min().copied()) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (a, b) => a.or(b),
        };
        for i in 0..list.len() {
            let before = i > 0 && gaps[i - 1] <= cluster_window;
            let after = i < gaps.len() && gaps[i] <= cluster_window;
            if before || after {
                s.clustered += 1;
            }
        }
    }
    Ok(s)
}


fn json_string(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn json_float(x: Option<f64>) -> String {
    match x {
        Some(v) if v.is_finite() => format!("{v:.6}"),
        _ => "null".to_string(),
    }
}

fn json_map<K: std::fmt::Display>(map: &BTreeMap<K, u64>) -> String {
    let items: Vec<String> = map.iter().map(|(k, v)| format!("{}: {v}", json_string(&k.to_string()))).collect();
    format!("{{{}}}", items.join(", "))
}

impl Stats {
    /// The statistics as one JSON object. A ratio that is undefined is `null`.
    pub fn to_json(&self) -> String {
        let contigs: Vec<String> = self
            .per_contig
            .iter()
            .map(|(name, c)| {
                format!(
                    "{}: {{\"length\": {}, \"count\": {}, \"per_kb\": {:.6}}}",
                    json_string(name),
                    c.length,
                    c.count,
                    c.per_kb
                )
            })
            .collect();
        let histogram: Vec<String> = self.fraction_histogram.iter().map(|n| n.to_string()).collect();
        let fields = [
            format!("\"total\": {}", self.total),
            format!("\"by_class\": {}", json_map(&self.by_class)),
            format!("\"by_effect\": {}", json_map(&self.by_effect)),
            format!("\"transitions\": {}", self.transitions),
            format!("\"transversions\": {}", self.transversions),
            format!("\"ts_tv\": {}", json_float(self.ts_tv)),
            format!("\"spectrum\": {}", json_map(&self.spectrum)),
            format!("\"spectrum_collapsed\": {}", json_map(&self.spectrum_collapsed)),
            format!("\"indel_lengths\": {}", json_map(&self.indel_lengths)),
            format!("\"frameshift_length_fraction\": {}", json_float(self.frameshift_length_fraction)),
            format!("\"per_contig\": {{{}}}", contigs.join(", ")),
            format!("\"overall_per_kb\": {}", json_float(self.overall_per_kb)),
            format!("\"fraction_histogram\": [{}]", histogram.join(", ")),
            format!("\"mid_fraction_calls\": {}", self.mid_fraction_calls),
            format!("\"mid_fraction_share\": {}", json_float(self.mid_fraction_share)),
            format!(
                "\"smallest_gap\": {}",
                self.smallest_gap.map(|g| g.to_string()).unwrap_or_else(|| "null".to_string())
            ),
            format!("\"clustered\": {}", self.clustered),
        ];
        format!("{{\n  {}\n}}\n", fields.join(",\n  "))
    }

    /// The statistics as a table with the columns `metric`, `key` and `value`.
    pub fn to_tsv(&self) -> String {
        let mut rows: Vec<(String, String, String)> = Vec::new();
        let scalar = |name: &str, value: String| (name.to_string(), String::new(), value);
        let ratio = |x: Option<f64>| x.map(|v| format!("{v:.4}")).unwrap_or_else(|| "NA".to_string());
        rows.push(scalar("total", self.total.to_string()));
        for (k, v) in &self.by_class {
            rows.push(("by_class".to_string(), k.clone(), v.to_string()));
        }
        for (k, v) in &self.by_effect {
            rows.push(("by_effect".to_string(), k.clone(), v.to_string()));
        }
        rows.push(scalar("transitions", self.transitions.to_string()));
        rows.push(scalar("transversions", self.transversions.to_string()));
        rows.push(scalar("ts_tv", ratio(self.ts_tv)));
        for (k, v) in &self.spectrum {
            rows.push(("spectrum".to_string(), k.clone(), v.to_string()));
        }
        for (k, v) in &self.spectrum_collapsed {
            rows.push(("spectrum_collapsed".to_string(), k.clone(), v.to_string()));
        }
        for (k, v) in &self.indel_lengths {
            rows.push(("indel_length".to_string(), k.to_string(), v.to_string()));
        }
        rows.push(scalar("frameshift_length_fraction", ratio(self.frameshift_length_fraction)));
        for (name, c) in &self.per_contig {
            rows.push(("per_contig_count".to_string(), name.clone(), c.count.to_string()));
            rows.push(("per_contig_per_kb".to_string(), name.clone(), format!("{:.4}", c.per_kb)));
        }
        rows.push(scalar("overall_per_kb", ratio(self.overall_per_kb)));
        for (i, n) in self.fraction_histogram.iter().enumerate() {
            let label = format!("{:.1}-{:.1}", i as f64 / 10.0, (i + 1) as f64 / 10.0);
            rows.push(("fraction_histogram".to_string(), label, n.to_string()));
        }
        rows.push(scalar("mid_fraction_calls", self.mid_fraction_calls.to_string()));
        rows.push(scalar("mid_fraction_share", ratio(self.mid_fraction_share)));
        rows.push(scalar("smallest_gap", self.smallest_gap.map(|g| g.to_string()).unwrap_or_else(|| "NA".to_string())));
        rows.push(scalar("clustered", self.clustered.to_string()));

        let mut out = String::from("metric\tkey\tvalue\n");
        for (metric, key, value) in rows {
            out.push_str(&format!("{metric}\t{key}\t{value}\n"));
        }
        out
    }
}

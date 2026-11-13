//! Combine the results of several strains (one VCF and one coverage BED each) into one table.
//!
//! A row is one variant: the same contig, position, REF and ALT. Each strain gets a cell for the row:
//! a call (with its allele fraction and depth), `weak` (no call, but the evidence file of the strain has some reads for this
//! ALT, below the calling thresholds), `ref` (no call, no such reads, and the coverage there is enough to say the base
//! matches the reference), `low` (no call, but too little coverage to say), `unknown` (no call and no coverage data), or
//! `other` (the strain has a different ALT at the same position).
//!
//! One strain can be marked as the parent. Then each other strain gets a status for the row: `new` (a call that the
//! parent lacks, and the parent is covered), `inherited` (the parent has it too), `lost` (the parent has it, the strain
//! does not, and the strain is covered), or `unknown` (the parent or the strain has weak evidence or too little coverage
//! to say). A variant is never `new` when the parent has some reads for it, and never `lost` when the strain has some.
//!
//! The evidence file of a strain is a VCF made with relaxed thresholds (`Caller::evidence_calls`). It is optional.
//!
//! The coverage BED has windows (`contig`, 0-based start, end, mean depth), as `depth-out` of the caller or mosdepth
//! writes it. A window is coarse: the depth at one site can differ from the window mean.

use std::collections::{BTreeMap, HashMap, HashSet};

/// The input of one strain. `depth_bed` is optional; without it, a missing call is `unknown`.
pub struct CombineInput<'a> {
    pub name: &'a str,
    pub vcf: &'a str,
    pub depth_bed: Option<&'a str>,
    /// Calls at relaxed thresholds, from the same mapping. Without it, a weak signal cannot be told from no signal.
    pub evidence_vcf: Option<&'a str>,
}

pub struct CombineOptions<'a> {
    /// The name of the parent strain, one of the inputs.
    pub parent: Option<&'a str>,
    /// A window with at least this mean depth counts as covered.
    pub min_depth: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Cell {
    Call { af: f64, depth: u64 },
    /// Some reads carry this ALT, but too few for a call.
    Weak { af: f64, depth: u64 },
    Ref,
    Low,
    Unknown,
    Other,
}

/// The status of a strain's cell compared with the parent's cell, for one row.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    /// This strain is the parent.
    Parent,
    New,
    Inherited,
    Lost,
    Unknown,
    /// Neither the strain nor the parent has the variant (another strain has it), or there is no parent.
    None,
}

#[derive(Debug)]
pub struct Row {
    pub contig: String,
    pub pos: u64,
    pub reference: String,
    pub alt: String,
    pub class: String,
    pub cells: Vec<Cell>,
    pub status: Vec<Status>,
}

#[derive(Debug)]
pub struct Combined {
    pub strains: Vec<String>,
    pub parent: Option<usize>,
    pub rows: Vec<Row>,
    /// Things the user should know that are not errors, such as a coverage contig that no VCF has.
    pub warnings: Vec<String>,
    /// The contigs with a known length, in the order of the first VCF that names them.
    pub contigs: Vec<(String, u64)>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct StrainSummary {
    pub strain: String,
    pub calls: usize,
    pub new: usize,
    pub inherited: usize,
    pub lost: usize,
    pub unknown: usize,
}

struct Record {
    contig: String,
    pos: u64,
    reference: String,
    alt: String,
    class: String,
    af: f64,
    depth: u64,
}

type ContigLengths = Vec<(String, u64)>;
type Windows = HashMap<String, Vec<(u64, u64, f64)>>;

fn parse_vcf(text: &str, name: &str) -> Result<(ContigLengths, Vec<Record>), String> {
    let mut contigs = Vec::new();
    let mut records = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("##contig=<") {
            let attr = |key: &str| rest.trim_end_matches('>').split(',').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='));
            if let (Some(id), Some(len)) = (attr("ID"), attr("length").and_then(|l| l.parse().ok())) {
                contigs.push((id.to_string(), len));
            }
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            return Err(format!("{name}: a VCF line has {} fields, needs at least 8: {line}", f.len()));
        }
        let pos: u64 = f[1].parse().map_err(|_| format!("{name}: bad VCF position in: {line}"))?;
        if pos == 0 {
            return Err(format!("{name}: VCF position 0 is not valid in: {line}"));
        }
        if !matches!(f[6], "PASS" | "." | "") {
            continue; // the record failed a filter
        }
        let info = |key: &str| f[7].split(';').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='));
        let depth: u64 = info("DP").and_then(|v| v.parse().ok()).unwrap_or(0);
        let fractions: Vec<&str> = info("AF").map(|v| v.split(',').collect()).unwrap_or_default();
        let counts: Vec<&str> = info("AO").map(|v| v.split(',').collect()).unwrap_or_default();
        for (k, alt) in f[4].split(',').enumerate() {
            if alt.is_empty() || alt == "*" {
                continue;
            }
            let af = fractions
                .get(k)
                .and_then(|v| v.parse::<f64>().ok())
                .or_else(|| (depth > 0).then(|| counts.get(k)?.parse::<f64>().ok().map(|n| n / depth as f64)).flatten())
                .unwrap_or(f64::NAN);
            let class = match info("CLASS") {
                Some(c) => c.to_string(),
                None => match crate::classify(f[3], alt, 50) {
                    crate::VariantClass::Snp => "SNP",
                    crate::VariantClass::SmallIndel => "SMALL_INDEL",
                    crate::VariantClass::LargeIndel => "LARGE_INDEL",
                }
                .to_string(),
            };
            records.push(Record { contig: f[0].to_string(), pos, reference: f[3].to_string(), alt: alt.to_string(), class, af, depth });
        }
    }
    Ok((contigs, records))
}

fn parse_bed(text: &str, name: &str) -> Result<Windows, String> {
    let mut windows: Windows = HashMap::new();
    for line in text.lines() {
        if line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            return Err(format!("{name}: a coverage line has {} fields, needs at least 4: {line}", f.len()));
        }
        let number = |text: &str| text.parse::<f64>().map_err(|_| format!("{name}: bad number in the coverage line: {line}"));
        windows.entry(f[0].to_string()).or_default().push((number(f[1])? as u64, number(f[2])? as u64, number(f[3])?));
    }
    for list in windows.values_mut() {
        list.sort_by_key(|w| w.0);
    }
    Ok(windows)
}

/// The mean depth of the window that holds the 1-based position, if there is one.
fn depth_at(windows: &Windows, contig: &str, pos: u64) -> Option<f64> {
    let list = windows.get(contig)?;
    let zero_based = pos - 1;
    let after = list.partition_point(|w| w.0 <= zero_based);
    let window = list.get(after.checked_sub(1)?)?;
    (zero_based < window.1).then_some(window.2)
}

fn status_of(cell: Cell, parent: Cell) -> Status {
    use Cell::*;
    match (cell, parent) {
        (Call { .. }, Call { .. }) => Status::Inherited,
        (Call { .. }, Ref | Other) => Status::New,
        (Call { .. }, Low | Unknown | Weak { .. }) => Status::Unknown,
        (Ref | Other, Call { .. }) => Status::Lost,
        (Low | Unknown | Weak { .. }, Call { .. }) => Status::Unknown,
        _ => Status::None,
    }
}

pub fn combine(inputs: &[CombineInput], options: &CombineOptions) -> Result<Combined, String> {
    if inputs.is_empty() {
        return Err("there are no strains to combine".to_string());
    }
    let mut seen = HashSet::new();
    for input in inputs {
        if input.name.is_empty() || input.name.contains(['\t', '\n', '\r']) {
            return Err(format!("the strain name {:?} is not valid: it must not be empty or hold a tab or a line break", input.name));
        }
        if !seen.insert(input.name) {
            return Err(format!("the strain name {} is used more than once", input.name));
        }
    }
    let parent = match options.parent {
        None => None,
        Some(name) => Some(
            inputs.iter().position(|i| i.name == name).ok_or_else(|| format!("the parent {name} is not one of the strains"))?,
        ),
    };

    let mut order: Vec<String> = Vec::new();
    let mut lengths: HashMap<String, (u64, &str)> = HashMap::new();
    let mut parsed = Vec::new();
    let mut windows = Vec::new();
    let mut evidence: Vec<HashMap<(String, u64, String, String), (f64, u64)>> = Vec::new();
    for input in inputs {
        let (contigs, records) = parse_vcf(input.vcf, input.name)?;
        for (contig, length) in contigs {
            match lengths.get(&contig) {
                Some((other, from)) if *other != length => {
                    return Err(format!("contig {contig} has length {other} in {from} and {length} in {}", input.name));
                }
                Some(_) => {}
                None => {
                    lengths.insert(contig.clone(), (length, input.name));
                    order.push(contig);
                }
            }
        }
        for r in &records {
            if !order.contains(&r.contig) {
                order.push(r.contig.clone());
            }
        }
        parsed.push(records);
        windows.push(input.depth_bed.map(|b| parse_bed(b, input.name)).transpose()?);
        let weak = match input.evidence_vcf {
            Some(text) => parse_vcf(text, input.name)?.1,
            None => Vec::new(),
        };
        evidence.push(weak.into_iter().map(|r| ((r.contig, r.pos, r.reference, r.alt), (r.af, r.depth))).collect());
    }
    let mut warnings = Vec::new();
    for (input, w) in inputs.iter().zip(&windows) {
        let mut unknown: Vec<&String> = w.iter().flat_map(|w| w.keys()).filter(|c| !order.contains(c)).collect();
        unknown.sort();
        for contig in unknown {
            warnings.push(format!("{}: the coverage file has contig {contig}, which no VCF has, so it is not used", input.name));
        }
    }
    let rank: HashMap<&str, usize> = order.iter().enumerate().map(|(i, c)| (c.as_str(), i)).collect();

    // Every distinct variant, with the call of each strain that has it.
    let mut variants: BTreeMap<(usize, u64, String, String), (String, Vec<Option<(f64, u64)>>)> = BTreeMap::new();
    let mut called_at: Vec<HashSet<(&str, u64)>> = vec![HashSet::new(); inputs.len()];
    for (i, records) in parsed.iter().enumerate() {
        for r in records {
            let key = (rank[r.contig.as_str()], r.pos, r.reference.clone(), r.alt.clone());
            let entry = variants.entry(key).or_insert_with(|| (r.class.clone(), vec![None; inputs.len()]));
            entry.1[i] = Some((r.af, r.depth));
            called_at[i].insert((r.contig.as_str(), r.pos));
        }
    }

    let mut rows = Vec::new();
    for ((contig_rank, pos, reference, alt), (class, calls)) in variants {
        let contig = order[contig_rank].clone();
        let cells: Vec<Cell> = (0..inputs.len())
            .map(|i| match calls[i] {
                Some((af, depth)) => Cell::Call { af, depth },
                None if called_at[i].contains(&(contig.as_str(), pos)) => Cell::Other,
                None if evidence[i].contains_key(&(contig.clone(), pos, reference.clone(), alt.clone())) => {
                    let (af, depth) = evidence[i][&(contig.clone(), pos, reference.clone(), alt.clone())];
                    Cell::Weak { af, depth }
                }
                None => match windows[i].as_ref().and_then(|w| depth_at(w, &contig, pos)) {
                    None => Cell::Unknown,
                    Some(d) if d >= options.min_depth => Cell::Ref,
                    Some(_) => Cell::Low,
                },
            })
            .collect();
        let status = (0..inputs.len())
            .map(|i| match parent {
                None => Status::None,
                Some(p) if p == i => Status::Parent,
                Some(p) => status_of(cells[i], cells[p]),
            })
            .collect();
        rows.push(Row { contig, pos, reference, alt, class, cells, status });
    }
    let contigs: Vec<(String, u64)> = order.iter().filter_map(|c| lengths.get(c).map(|(len, _)| (c.clone(), *len))).collect();
    Ok(Combined { strains: inputs.iter().map(|i| i.name.to_string()).collect(), parent, rows, warnings, contigs })
}

impl Cell {
    fn text(&self) -> String {
        match self {
            Cell::Call { af, depth } => format!("{af:.2}:{depth}"),
            Cell::Weak { af, depth } => format!("weak:{af:.2}:{depth}"),
            Cell::Ref => "ref".to_string(),
            Cell::Low => "low".to_string(),
            Cell::Unknown => "unknown".to_string(),
            Cell::Other => "other".to_string(),
        }
    }
}

impl Status {
    fn text(&self) -> &'static str {
        match self {
            Status::Parent => "parent",
            Status::New => "new",
            Status::Inherited => "inherited",
            Status::Lost => "lost",
            Status::Unknown => "unknown",
            Status::None => ".",
        }
    }
}

/// A sample name for an INFO value: `,` `;` `=` space `%` and control characters become %XX, so the names can be listed with commas.
fn escape_info(name: &str) -> String {
    // Work on bytes: only ASCII bytes are changed, so the UTF-8 sequences of other characters stay intact.
    let mut bytes = Vec::with_capacity(name.len());
    for b in name.bytes() {
        if matches!(b, b',' | b';' | b'=' | b' ' | b'%') || b < 0x20 || b == 0x7f {
            bytes.extend_from_slice(format!("%{b:02X}").as_bytes());
        } else {
            bytes.push(b);
        }
    }
    String::from_utf8(bytes).expect("escaping ASCII bytes keeps valid UTF-8")
}

impl Combined {
    /// One VCF record per distinct variant: the samples that call it (`NS`, `SAMPLES`), the samples where it is new compared with the
    /// parent (`NEW`), and the mean AF and depth of the calls. The record can go through the summary, the annotation and the locus
    /// drawing like the VCF of one sample.
    pub fn to_vcf(&self) -> String {
        let mut out = String::from("##fileformat=VCFv4.2\n");
        for (name, len) in &self.contigs {
            out.push_str(&format!("##contig=<ID={name},length={len}>\n"));
        }
        out.push_str("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with a call\">\n");
        out.push_str("##INFO=<ID=SAMPLES,Number=.,Type=String,Description=\"Samples with a call (a name has , ; = space and % written as %XX)\">\n");
        out.push_str("##INFO=<ID=NEW,Number=.,Type=String,Description=\"Samples where the variant is new compared with the parent\">\n");
        out.push_str("##INFO=<ID=CLASS,Number=1,Type=String,Description=\"SNP, SMALL_INDEL or LARGE_INDEL\">\n");
        out.push_str("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Mean depth of the calls\">\n");
        out.push_str("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Mean alternate allele fraction of the calls\">\n");
        out.push_str(&format!("##qmutfind_samples={}\n", self.strains.len()));
        if let Some(p) = self.parent {
            out.push_str(&format!("##qmutfind_parent={}\n", escape_info(&self.strains[p])));
        }
        out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for r in &self.rows {
            let calls: Vec<(usize, f64, u64)> = r
                .cells
                .iter()
                .enumerate()
                .filter_map(|(i, c)| if let Cell::Call { af, depth } = c { Some((i, *af, *depth)) } else { None })
                .collect();
            let names = calls.iter().map(|(i, _, _)| escape_info(&self.strains[*i])).collect::<Vec<_>>().join(",");
            let mut info = format!("NS={};SAMPLES={names};CLASS={}", calls.len(), r.class);
            let depth = calls.iter().map(|c| c.2 as f64).sum::<f64>() / calls.len().max(1) as f64;
            info.push_str(&format!(";DP={}", depth.round() as u64));
            let afs: Vec<f64> = calls.iter().map(|c| c.1).filter(|a| a.is_finite()).collect();
            if !afs.is_empty() {
                info.push_str(&format!(";AF={:.4}", afs.iter().sum::<f64>() / afs.len() as f64));
            }
            let new: Vec<String> = r.status.iter().enumerate().filter(|(_, s)| **s == Status::New).map(|(i, _)| escape_info(&self.strains[i])).collect();
            if !new.is_empty() {
                info.push_str(&format!(";NEW={}", new.join(",")));
            }
            out.push_str(&format!("{}\t{}\t.\t{}\t{}\t.\tPASS\t{info}\n", r.contig, r.pos, r.reference, r.alt));
        }
        out
    }

    /// One line per variant. The cells come first, then (with a parent) one status column per strain.
    pub fn to_tsv(&self) -> String {
        let mut header = vec!["CHROM".to_string(), "POS".into(), "REF".into(), "ALT".into(), "CLASS".into(), "strains".into()];
        header.extend(self.strains.iter().cloned());
        if self.parent.is_some() {
            header.extend(self.strains.iter().map(|s| format!("{s}.status")));
        }
        let mut out = header.join("\t") + "\n";
        for r in &self.rows {
            let called = r.cells.iter().filter(|c| matches!(c, Cell::Call { .. })).count();
            let mut fields = vec![r.contig.clone(), r.pos.to_string(), r.reference.clone(), r.alt.clone(), r.class.clone(), called.to_string()];
            fields.extend(r.cells.iter().map(Cell::text));
            if self.parent.is_some() {
                fields.extend(r.status.iter().map(|s| s.text().to_string()));
            }
            out.push_str(&fields.join("\t"));
            out.push('\n');
        }
        out
    }

    /// Counts per strain. `new`, `inherited`, `lost` and `unknown` are zero without a parent.
    pub fn summary(&self) -> Vec<StrainSummary> {
        let mut out: Vec<StrainSummary> = self
            .strains
            .iter()
            .map(|s| StrainSummary { strain: s.clone(), calls: 0, new: 0, inherited: 0, lost: 0, unknown: 0 })
            .collect();
        for r in &self.rows {
            for (i, summary) in out.iter_mut().enumerate() {
                if matches!(r.cells[i], Cell::Call { .. }) {
                    summary.calls += 1;
                }
                match r.status[i] {
                    Status::New => summary.new += 1,
                    Status::Inherited => summary.inherited += 1,
                    Status::Lost => summary.lost += 1,
                    Status::Unknown => summary.unknown += 1,
                    Status::Parent | Status::None => {}
                }
            }
        }
        out
    }

    pub fn summary_tsv(&self) -> String {
        let mut out = String::from("strain\tcalls\tnew\tinherited\tlost\tunknown\n");
        for s in self.summary() {
            out.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\n", s.strain, s.calls, s.new, s.inherited, s.lost, s.unknown));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vcf(contigs: &[(&str, u64)], records: &[(&str, u64, &str, &str, &str, u64, f64)]) -> String {
        let mut text = String::from("##fileformat=VCFv4.2\n");
        for (name, len) in contigs {
            text.push_str(&format!("##contig=<ID={name},length={len}>\n"));
        }
        text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for (contig, pos, reference, alt, class, dp, af) in records {
            let ao = (*dp as f64 * af).round() as u64;
            text.push_str(&format!("{contig}\t{pos}\t.\t{reference}\t{alt}\t50.0\tPASS\tSAMPLE=s;CLASS={class};DP={dp};AO={ao};AF={af:.4}\n"));
        }
        text
    }

    /// Windows of 1,000 bases on `chrA` (length 3,000) with the given mean depths.
    fn bed(depths: [f64; 3]) -> String {
        depths.iter().enumerate().map(|(i, d)| format!("chrA\t{}\t{}\t{d}\n", i * 1000, (i + 1) * 1000)).collect()
    }

    const CONTIGS: &[(&str, u64)] = &[("chrA", 3000)];

    fn options(parent: Option<&str>) -> CombineOptions<'_> {
        CombineOptions { parent, min_depth: 10.0 }
    }

    fn two_strains() -> (String, String, String) {
        // The parent has a SNP at 500. The evolved strain has that SNP, a new SNP at 1500 and a new indel at 2500.
        let parent = vcf(CONTIGS, &[("chrA", 500, "A", "G", "SNP", 30, 1.0)]);
        let evolved = vcf(CONTIGS, &[
            ("chrA", 500, "A", "G", "SNP", 28, 1.0),
            ("chrA", 1500, "C", "T", "SNP", 25, 0.96),
            ("chrA", 2500, "AT", "A", "SMALL_INDEL", 20, 0.9),
        ]);
        (parent, evolved, bed([30.0, 30.0, 30.0]))
    }

    fn run(parent_name: Option<&str>) -> Combined {
        let (p, e, b) = two_strains();
        let inputs = [
            CombineInput { name: "parent", vcf: &p, depth_bed: Some(&b), evidence_vcf: None },
            CombineInput { name: "evolved", vcf: &e, depth_bed: Some(&b), evidence_vcf: None },
        ];
        combine(&inputs, &options(parent_name)).unwrap()
    }

    #[test]
    fn one_row_per_distinct_variant_sorted_by_contig_order_then_position() {
        let c = run(None);

        let keys: Vec<(u64, &str, &str)> = c.rows.iter().map(|r| (r.pos, r.reference.as_str(), r.alt.as_str())).collect();
        assert_eq!(keys, vec![(500, "A", "G"), (1500, "C", "T"), (2500, "AT", "A")]);
        assert_eq!(c.rows[2].class, "SMALL_INDEL");
        assert_eq!(c.strains, vec!["parent", "evolved"]);
    }

    #[test]
    fn a_strain_with_the_call_gets_its_fraction_and_depth_and_one_without_it_is_ref_when_covered() {
        let c = run(None);

        assert_eq!(c.rows[1].cells[0], Cell::Ref); // the parent has no call at 1500, and window depth is 30
        assert_eq!(c.rows[1].cells[1], Cell::Call { af: 0.96, depth: 25 });
        assert_eq!(c.rows[0].cells[0], Cell::Call { af: 1.0, depth: 30 });
    }

    #[test]
    fn a_missing_call_is_low_when_the_window_is_below_the_minimum_and_unknown_without_a_bed() {
        let (p, e, _) = two_strains();
        let low = bed([30.0, 4.0, 30.0]); // the window of position 1500 has depth 4
        let inputs = [
            CombineInput { name: "parent", vcf: &p, depth_bed: Some(&low), evidence_vcf: None },
            CombineInput { name: "evolved", vcf: &e, depth_bed: None, evidence_vcf: None },
        ];

        let c = combine(&inputs, &options(None)).unwrap();

        assert_eq!(c.rows[1].cells[0], Cell::Low);
        assert_eq!(c.rows[0].cells[1], Cell::Call { af: 1.0, depth: 28 });
        let no_bed = combine(&[
            CombineInput { name: "parent", vcf: &p, depth_bed: None, evidence_vcf: None },
            CombineInput { name: "evolved", vcf: &e, depth_bed: None, evidence_vcf: None },
        ], &options(None)).unwrap();
        assert_eq!(no_bed.rows[1].cells[0], Cell::Unknown);
    }

    #[test]
    fn a_window_with_exactly_the_minimum_depth_counts_as_covered() {
        let (p, e, _) = two_strains();
        let edge = bed([30.0, 10.0, 30.0]); // 10.0 is the minimum depth of the options
        let inputs = [
            CombineInput { name: "parent", vcf: &p, depth_bed: Some(&edge), evidence_vcf: None },
            CombineInput { name: "evolved", vcf: &e, depth_bed: Some(&edge), evidence_vcf: None },
        ];

        let c = combine(&inputs, &options(None)).unwrap();

        assert_eq!(c.rows[1].cells[0], Cell::Ref);
    }

    #[test]
    fn a_different_alt_at_the_same_position_is_a_separate_row_and_the_other_strain_is_marked_other() {
        let a = vcf(CONTIGS, &[("chrA", 700, "A", "G", "SNP", 30, 1.0)]);
        let b = vcf(CONTIGS, &[("chrA", 700, "A", "T", "SNP", 30, 1.0)]);
        let bd = bed([30.0, 30.0, 30.0]);
        let inputs = [
            CombineInput { name: "a", vcf: &a, depth_bed: Some(&bd), evidence_vcf: None },
            CombineInput { name: "b", vcf: &b, depth_bed: Some(&bd), evidence_vcf: None },
        ];

        let c = combine(&inputs, &options(None)).unwrap();

        assert_eq!(c.rows.len(), 2);
        let g = c.rows.iter().find(|r| r.alt == "G").unwrap();
        assert_eq!(g.cells[1], Cell::Other);
        assert!(matches!(g.cells[0], Cell::Call { .. }));
    }

    #[test]
    fn with_a_parent_each_row_is_new_inherited_or_lost_for_every_other_strain() {
        let c = run(Some("parent"));

        assert_eq!(c.parent, Some(0));
        let status = |row: usize| c.rows[row].status.clone();
        assert_eq!(status(0), vec![Status::Parent, Status::Inherited]);
        assert_eq!(status(1), vec![Status::Parent, Status::New]);
        assert_eq!(status(2), vec![Status::Parent, Status::New]);

        // A strain that lost the parent's variant, and covered there.
        let (p, _, b) = two_strains();
        let reverted = vcf(CONTIGS, &[]);
        let inputs = [
            CombineInput { name: "parent", vcf: &p, depth_bed: Some(&b), evidence_vcf: None },
            CombineInput { name: "revertant", vcf: &reverted, depth_bed: Some(&b), evidence_vcf: None },
        ];
        let c = combine(&inputs, &options(Some("parent"))).unwrap();
        assert_eq!(c.rows[0].status, vec![Status::Parent, Status::Lost]);
    }

    #[test]
    fn a_parent_with_too_little_coverage_makes_the_status_unknown_instead_of_new() {
        let (p, e, _) = two_strains();
        let parent_bed = bed([30.0, 3.0, 30.0]); // the parent has depth 3 at 1500
        let strain_bed = bed([30.0, 30.0, 30.0]);
        let inputs = [
            CombineInput { name: "parent", vcf: &p, depth_bed: Some(&parent_bed), evidence_vcf: None },
            CombineInput { name: "evolved", vcf: &e, depth_bed: Some(&strain_bed), evidence_vcf: None },
        ];

        let c = combine(&inputs, &options(Some("parent"))).unwrap();

        assert_eq!(c.rows[1].status[1], Status::Unknown);
        assert_eq!(c.rows[2].status[1], Status::New);
    }

    #[test]
    fn with_three_strains_a_variant_neither_the_strain_nor_the_parent_has_has_no_status() {
        let (p, e, b) = two_strains();
        let other = vcf(CONTIGS, &[]);
        let inputs = [
            CombineInput { name: "parent", vcf: &p, depth_bed: Some(&b), evidence_vcf: None },
            CombineInput { name: "evolved", vcf: &e, depth_bed: Some(&b), evidence_vcf: None },
            CombineInput { name: "plain", vcf: &other, depth_bed: Some(&b), evidence_vcf: None },
        ];

        let c = combine(&inputs, &options(Some("parent"))).unwrap();

        assert_eq!(c.rows[1].status, vec![Status::Parent, Status::New, Status::None]);
    }

    #[test]
    fn the_input_is_checked_for_a_repeated_name_a_missing_parent_and_different_contig_lengths() {
        let (p, e, b) = two_strains();
        let inputs = [
            CombineInput { name: "x", vcf: &p, depth_bed: Some(&b), evidence_vcf: None },
            CombineInput { name: "x", vcf: &e, depth_bed: Some(&b), evidence_vcf: None },
        ];
        assert!(combine(&inputs, &options(None)).unwrap_err().contains("more than once"));

        let inputs = [
            CombineInput { name: "a", vcf: &p, depth_bed: Some(&b), evidence_vcf: None },
            CombineInput { name: "b", vcf: &e, depth_bed: Some(&b), evidence_vcf: None },
        ];
        assert!(combine(&inputs, &options(Some("zzz"))).unwrap_err().contains("parent"));

        let other_genome = vcf(&[("chrA", 9999)], &[]);
        let inputs = [
            CombineInput { name: "a", vcf: &p, depth_bed: None, evidence_vcf: None },
            CombineInput { name: "b", vcf: &other_genome, depth_bed: None, evidence_vcf: None },
        ];
        assert!(combine(&inputs, &options(None)).unwrap_err().contains("chrA"));

        assert!(combine(&[], &options(None)).is_err());
    }

    #[test]
    fn the_tsv_has_a_header_the_cells_and_with_a_parent_the_status_columns() {
        let c = run(Some("parent"));

        let tsv = c.to_tsv();
        let lines: Vec<&str> = tsv.lines().collect();

        assert_eq!(lines[0], "CHROM\tPOS\tREF\tALT\tCLASS\tstrains\tparent\tevolved\tparent.status\tevolved.status");
        assert_eq!(lines[1], "chrA\t500\tA\tG\tSNP\t2\t1.00:30\t1.00:28\tparent\tinherited");
        assert_eq!(lines[2], "chrA\t1500\tC\tT\tSNP\t1\tref\t0.96:25\tparent\tnew");
        assert_eq!(lines.len(), 4);
        let no_parent = run(None).to_tsv();
        assert_eq!(no_parent.lines().next().unwrap(), "CHROM\tPOS\tREF\tALT\tCLASS\tstrains\tparent\tevolved");
    }

    #[test]
    fn the_summary_counts_calls_and_the_status_of_each_strain() {
        let c = run(Some("parent"));

        let s = c.summary();

        assert_eq!(s[0], StrainSummary { strain: "parent".into(), calls: 1, new: 0, inherited: 0, lost: 0, unknown: 0 });
        assert_eq!(s[1], StrainSummary { strain: "evolved".into(), calls: 3, new: 2, inherited: 1, lost: 0, unknown: 0 });
        assert_eq!(c.summary_tsv().lines().next().unwrap(), "strain\tcalls\tnew\tinherited\tlost\tunknown");
    }

    fn raw_vcf(lines: &[&str]) -> String {
        let mut text = String::from("##fileformat=VCFv4.2\n##contig=<ID=chrA,length=3000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for l in lines {
            text.push_str(l);
            text.push('\n');
        }
        text
    }

    #[test]
    fn a_record_with_a_filter_other_than_pass_or_dot_is_ignored() {
        let v = raw_vcf(&[
            "chrA\t500\t.\tA\tG\t50\tPASS\tDP=30;AO=30;AF=1.0",
            "chrA\t600\t.\tA\tG\t50\t.\tDP=30;AO=30;AF=1.0",
            "chrA\t700\t.\tA\tG\t5\tLowQual\tDP=30;AO=3;AF=0.1",
        ]);

        let c = combine(&[CombineInput { name: "a", vcf: &v, depth_bed: None, evidence_vcf: None }], &options(None)).unwrap();

        let positions: Vec<u64> = c.rows.iter().map(|r| r.pos).collect();
        assert_eq!(positions, vec![500, 600]);
    }

    #[test]
    fn each_alt_of_a_multi_allelic_record_is_its_own_row_with_its_own_fraction() {
        let v = raw_vcf(&["chrA\t500\t.\tA\tG,T\t50\tPASS\tDP=30;AO=18,9;AF=0.6000,0.3000"]);

        let c = combine(&[CombineInput { name: "a", vcf: &v, depth_bed: None, evidence_vcf: None }], &options(None)).unwrap();

        assert_eq!(c.rows.len(), 2);
        assert_eq!((c.rows[0].alt.as_str(), c.rows[0].cells[0]), ("G", Cell::Call { af: 0.6, depth: 30 }));
        assert_eq!((c.rows[1].alt.as_str(), c.rows[1].cells[0]), ("T", Cell::Call { af: 0.3, depth: 30 }));
    }

    #[test]
    fn position_zero_a_bad_name_and_a_bed_for_an_unknown_contig_are_reported() {
        let zero = raw_vcf(&["chrA\t0\t.\tA\tG\t50\tPASS\tDP=30;AO=30;AF=1.0"]);
        let err = combine(&[CombineInput { name: "a", vcf: &zero, depth_bed: None, evidence_vcf: None }], &options(None)).unwrap_err();
        assert!(err.contains("position"), "{err}");

        let ok = raw_vcf(&[]);
        for bad in ["", "a\tb", "a\nb"] {
            let err = combine(&[CombineInput { name: bad, vcf: &ok, depth_bed: None, evidence_vcf: None }], &options(None)).unwrap_err();
            assert!(err.contains("name"), "{bad:?}: {err}");
        }

        let wrong_bed = "chrZ\t0\t1000\t30\n";
        let c = combine(&[CombineInput { name: "a", vcf: &ok, depth_bed: Some(wrong_bed), evidence_vcf: None }], &options(None)).unwrap();
        assert!(c.warnings.iter().any(|w| w.contains("chrZ") && w.contains("a")), "{:?}", c.warnings);
    }

    #[test]
    fn evidence_below_the_calling_thresholds_makes_a_weak_cell_that_blocks_new_and_lost() {
        let parent = vcf(CONTIGS, &[]);
        let parent_evidence = vcf(CONTIGS, &[("chrA", 1500, "C", "T", "SNP", 30, 0.2)]); // 6 of 30 reads, no call
        let child = vcf(CONTIGS, &[("chrA", 1500, "C", "T", "SNP", 25, 0.96)]);
        let b = bed([30.0, 30.0, 30.0]);
        let inputs = [
            CombineInput { name: "parent", vcf: &parent, depth_bed: Some(&b), evidence_vcf: Some(&parent_evidence) },
            CombineInput { name: "child", vcf: &child, depth_bed: Some(&b), evidence_vcf: None },
        ];

        let c = combine(&inputs, &options(Some("parent"))).unwrap();

        assert_eq!(c.rows[0].cells[0], Cell::Weak { af: 0.2, depth: 30 });
        assert_eq!(c.rows[0].status[1], Status::Unknown); // not "new": the parent has 20% alternate reads
        assert!(c.to_tsv().contains("weak:0.20:30"));
        // Without the evidence file the same site is "new".
        let inputs = [
            CombineInput { name: "parent", vcf: &parent, depth_bed: Some(&b), evidence_vcf: None },
            CombineInput { name: "child", vcf: &child, depth_bed: Some(&b), evidence_vcf: None },
        ];
        assert_eq!(combine(&inputs, &options(Some("parent"))).unwrap().rows[0].status[1], Status::New);
        // A weak cell in the child blocks "lost" where the parent has the call.
        let parent_call = vcf(CONTIGS, &[("chrA", 1500, "C", "T", "SNP", 25, 0.96)]);
        let child_evidence = vcf(CONTIGS, &[("chrA", 1500, "C", "T", "SNP", 30, 0.2)]);
        let inputs = [
            CombineInput { name: "parent", vcf: &parent_call, depth_bed: Some(&b), evidence_vcf: None },
            CombineInput { name: "child", vcf: &parent, depth_bed: Some(&b), evidence_vcf: Some(&child_evidence) },
        ];
        assert_eq!(combine(&inputs, &options(Some("parent"))).unwrap().rows[0].status[1], Status::Unknown);
    }

    #[test]
    fn the_merged_vcf_has_one_record_per_variant_with_the_number_of_samples_and_their_names() {
        let c = run(Some("parent"));

        let vcf = c.to_vcf();
        let lines: Vec<&str> = vcf.lines().collect();

        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        assert!(lines.contains(&"##contig=<ID=chrA,length=3000>"));
        assert!(lines.iter().any(|l| l.starts_with("##INFO=<ID=NS,")) && lines.iter().any(|l| l.starts_with("##INFO=<ID=SAMPLES,")));
        assert!(lines.contains(&"##qmutfind_samples=2"));
        let records: Vec<Vec<&str>> = lines.iter().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').collect()).collect();
        assert_eq!(records.len(), 3);
        // Position 500 is called in both strains: its AF is the mean (1.00), DP the mean (29), and the status of the child is inherited.
        assert_eq!(&records[0][..5], &["chrA", "500", ".", "A", "G"]);
        assert!(records[0][7].contains("NS=2") && records[0][7].contains("SAMPLES=parent,evolved") && records[0][7].contains("CLASS=SNP"), "{}", records[0][7]);
        assert!(records[0][7].contains("DP=29") && records[0][7].contains("AF=1.0000"), "{}", records[0][7]);
        assert!(!records[0][7].contains("NEW="), "{}", records[0][7]);
        // Position 1500 is called in one strain and is new there.
        assert!(records[1][7].contains("NS=1") && records[1][7].contains("SAMPLES=evolved") && records[1][7].contains("NEW=evolved"), "{}", records[1][7]);
    }

    #[test]
    fn the_merged_vcf_without_a_parent_has_no_new_tag_and_names_are_escaped() {
        let v = vcf(CONTIGS, &[("chrA", 500, "A", "G", "SNP", 30, 1.0)]);
        let b = bed([30.0, 30.0, 30.0]);
        let inputs = [
            CombineInput { name: "a,b;c=d %e", vcf: &v, depth_bed: Some(&b), evidence_vcf: None },
            CombineInput { name: "plain", vcf: &v, depth_bed: Some(&b), evidence_vcf: None },
        ];

        let c = combine(&inputs, &options(None)).unwrap();
        let vcf = c.to_vcf();

        let record = vcf.lines().find(|l| !l.starts_with('#')).unwrap();
        assert!(record.contains("SAMPLES=a%2Cb%3Bc%3Dd%20%25e,plain"), "{record}");
        assert!(!record.contains("NEW="));
        assert!(record.contains("NS=2"));
    }

    #[test]
    fn the_merged_vcf_keeps_a_variant_called_in_one_strain_even_when_the_others_have_weak_evidence() {
        let parent = vcf(CONTIGS, &[]);
        let weak = vcf(CONTIGS, &[("chrA", 1500, "C", "T", "SNP", 30, 0.2)]);
        let child = vcf(CONTIGS, &[("chrA", 1500, "C", "T", "SNP", 25, 0.96)]);
        let b = bed([30.0, 30.0, 30.0]);
        let inputs = [
            CombineInput { name: "parent", vcf: &parent, depth_bed: Some(&b), evidence_vcf: Some(&weak) },
            CombineInput { name: "child", vcf: &child, depth_bed: Some(&b), evidence_vcf: None },
        ];

        let vcf = combine(&inputs, &options(Some("parent"))).unwrap().to_vcf();

        let record = vcf.lines().find(|l| !l.starts_with('#')).unwrap();
        assert!(record.contains("NS=1") && record.contains("SAMPLES=child"), "{record}");
        assert!(!record.contains("NEW="), "the parent has weak evidence, so the status is unknown and not new: {record}");
    }

    #[test]
    fn a_name_with_non_ascii_letters_stays_readable_in_the_merged_vcf() {
        let v = vcf(CONTIGS, &[("chrA", 500, "A", "G", "SNP", 30, 1.0)]);
        let inputs = [CombineInput { name: "souche é1", vcf: &v, depth_bed: None, evidence_vcf: None }];

        let vcf_text = combine(&inputs, &options(None)).unwrap().to_vcf();

        assert!(vcf_text.contains("SAMPLES=souche%20é1"), "{vcf_text}");
    }
}

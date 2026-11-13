//! Effect of a variant on the canonical transcript of a gene.
//!
//! The GFF3 is read into transcripts. UTRs and introns are derived from the exons and the CDS
//! when the file has no such features. The canonical transcript of a gene is the one with the
//! longest CDS. Positions are 1-based and inclusive, as in GFF and VCF.

use std::collections::HashMap;
use std::io::{self, BufRead};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

#[derive(Debug, Clone)]
pub struct Transcript {
    pub id: String,
    pub gene: String,
    pub contig: String,
    pub strand: Strand,
    /// Exons, sorted by position.
    pub exons: Vec<(u64, u64)>,
    /// CDS blocks, sorted by position.
    pub cds: Vec<(u64, u64)>,
    pub transl_table: u8,
    /// The gene's name (for example `hrmA`) when the GFF gives one that is not just the locus tag.
    pub gene_name: Option<String>,
    /// The product (for example `response regulator`), from the mRNA or else the CDS line.
    pub product: Option<String>,
    /// The nucleotide accession of the transcript (for example `XM_001481640.2`) when the mRNA line of the GFF gives one. See `nucleotide_accession`.
    pub accession: Option<String>,
}

impl Transcript {
    fn start(&self) -> u64 {
        self.exons.iter().map(|e| e.0).min().unwrap_or(0)
    }

    fn end(&self) -> u64 {
        self.exons.iter().map(|e| e.1).max().unwrap_or(0)
    }

    fn cds_len(&self) -> u64 {
        self.cds.iter().map(|c| c.1 - c.0 + 1).sum()
    }

    /// Index of `pos` along the coding sequence (0-based, in transcript direction), if it is in the CDS.
    pub(crate) fn coding_index(&self, pos: u64) -> Option<u64> {
        let mut before = 0;
        match self.strand {
            Strand::Plus => {
                for &(s, e) in &self.cds {
                    if pos >= s && pos <= e {
                        return Some(before + pos - s);
                    }
                    before += e - s + 1;
                }
            }
            Strand::Minus => {
                for &(s, e) in self.cds.iter().rev() {
                    if pos >= s && pos <= e {
                        return Some(before + e - pos);
                    }
                    before += e - s + 1;
                }
            }
        }
        None
    }

    /// The coding sequence in transcript direction (reverse complemented on the minus strand).
    pub(crate) fn coding_sequence(&self, genome: &Genome) -> Option<Vec<u8>> {
        let mut out = Vec::with_capacity(self.cds_len() as usize);
        for idx in 0..self.cds_len() {
            let base = genome.base(&self.contig, self.position_of(idx)?)?.to_ascii_uppercase();
            out.push(if self.strand == Strand::Minus { complement(base) } else { base });
        }
        Some(out)
    }

    /// Genomic position of coding index `idx`.
    fn position_of(&self, idx: u64) -> Option<u64> {
        let mut left = idx;
        match self.strand {
            Strand::Plus => {
                for &(s, e) in &self.cds {
                    let len = e - s + 1;
                    if left < len {
                        return Some(s + left);
                    }
                    left -= len;
                }
            }
            Strand::Minus => {
                for &(s, e) in self.cds.iter().rev() {
                    let len = e - s + 1;
                    if left < len {
                        return Some(e - left);
                    }
                    left -= len;
                }
            }
        }
        None
    }
}

/// All canonical transcripts, by contig, and the number of transcripts of each gene.
pub struct Annotation {
    by_contig: HashMap<String, Vec<Transcript>>,
    transcripts_per_gene: HashMap<String, usize>,
}

fn attribute<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    attrs.split(';').find_map(|kv| kv.trim().strip_prefix(key)?.strip_prefix('='))
}

/// Whether `text` has the shape of an NCBI nucleotide accession with an optional version: one or two capital letters, an optional
/// underscore, 5 to 12 digits, then optionally `.` and digits (`XM_001481640.2`, `NM_000001.3`, `AB123456.1`).
fn is_accession(text: &str) -> bool {
    let (body, version) = match text.split_once('.') {
        Some((body, version)) => (body, Some(version)),
        None => (text, None),
    };
    if let Some(v) = version {
        if v.is_empty() || !v.bytes().all(|b| b.is_ascii_digit()) {
            return false;
        }
    }
    let letters = body.bytes().take_while(|b| b.is_ascii_uppercase()).count();
    let rest = &body[letters..];
    let rest = rest.strip_prefix('_').unwrap_or(rest);
    (1..=2).contains(&letters) && (5..=12).contains(&rest.len()) && rest.bytes().all(|b| b.is_ascii_digit())
}

/// The NCBI nucleotide accession of the transcript on a GFF3 mRNA line: a `GenBank` entry in `Dbxref` (any capitals; RefSeq files write
/// `GenBank:XM_...`), or else a `transcript_id` that is an accession. Nothing is guessed: a GFF without either gives `None`.
fn nucleotide_accession(attrs: &str) -> Option<String> {
    attribute(attrs, "Dbxref")
        .and_then(|d| {
            d.split(',').find_map(|entry| {
                let (source, id) = entry.split_once(':')?;
                (source.eq_ignore_ascii_case("genbank") && is_accession(id)).then(|| id.to_string())
            })
        })
        .or_else(|| attribute(attrs, "transcript_id").filter(|t| is_accession(t)).map(str::to_string))
}

/// Decode the `%XX` escapes of a GFF3 attribute value. A `%` that is not followed by two hex digits stays as it is.
fn percent_decode(text: &str) -> String {
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

/// Encode a value for a VCF INFO field: the characters that would end or break the field become `%XX`.
pub(crate) fn percent_encode_info(text: &str) -> String {
    let mut out = String::new();
    for c in text.chars() {
        match c {
            '%' | ';' | '=' | ',' | ' ' | '\t' | '\n' | '\r' => out.push_str(&format!("%{:02X}", c as u32)),
            other => out.push(other),
        }
    }
    out
}

struct Feature {
    kind: String,
    contig: String,
    strand: Strand,
    parent: Option<String>,
    label: String,
    name: Option<String>,
    product: Option<String>,
    accession: Option<String>,
}

impl Annotation {
    /// An annotation with no genes.
    pub fn empty() -> Annotation {
        Annotation { by_contig: HashMap::new(), transcripts_per_gene: HashMap::new() }
    }

    /// The canonical transcripts that overlap `contig:start-end`.
    pub fn transcripts_in(&self, contig: &str, start: u64, end: u64) -> Vec<&Transcript> {
        self.by_contig
            .get(contig)
            .map(|txs| txs.iter().filter(|t| t.start() <= end && t.end() >= start).collect())
            .unwrap_or_default()
    }

    /// Read a GFF3 file.
    pub fn from_gff<R: BufRead>(reader: R) -> io::Result<Annotation> {
        let mut features: HashMap<String, Feature> = HashMap::new();
        let mut exons: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let mut cds: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let mut tables: HashMap<String, u8> = HashMap::new();
        let mut cds_products: HashMap<String, String> = HashMap::new();
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 9 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, format!("GFF line has {} fields: {line}", f.len())));
            }
            let bad = |what: &str| io::Error::new(io::ErrorKind::InvalidData, format!("bad {what} in GFF line: {line}"));
            let start: u64 = f[3].parse().map_err(|_| bad("start"))?;
            let end: u64 = f[4].parse().map_err(|_| bad("end"))?;
            let strand = if f[6] == "-" { Strand::Minus } else { Strand::Plus };
            let attrs = f[8];
            match f[2] {
                "exon" | "CDS" => {
                    for parent in attribute(attrs, "Parent").unwrap_or("").split(',').filter(|p| !p.is_empty()) {
                        if f[2] == "exon" {
                            exons.entry(parent.to_string()).or_default().push((start, end));
                        } else {
                            cds.entry(parent.to_string()).or_default().push((start, end));
                            if let Some(t) = attribute(attrs, "transl_table").and_then(|t| t.parse().ok()) {
                                tables.insert(parent.to_string(), t);
                            }
                            if let Some(p) = attribute(attrs, "product") {
                                cds_products.entry(parent.to_string()).or_insert_with(|| percent_decode(p));
                            }
                        }
                    }
                }
                kind => {
                    let Some(id) = attribute(attrs, "ID") else { continue };
                    let label = attribute(attrs, "locus_tag")
                        .or_else(|| attribute(attrs, "Name"))
                        .unwrap_or(id)
                        .to_string();
                    // A gene's name: the `gene` attribute, or for a gene line a `Name` that is not the locus tag itself.
                    let name = attribute(attrs, "gene")
                        .map(percent_decode)
                        .or_else(|| if matches!(kind, "gene" | "pseudogene") { attribute(attrs, "Name").map(percent_decode) } else { None })
                        .filter(|n| !n.is_empty() && *n != label);
                    features.insert(
                        id.to_string(),
                        Feature {
                            kind: kind.to_string(),
                            contig: f[0].to_string(),
                            strand,
                            parent: attribute(attrs, "Parent").map(|p| p.split(',').next().unwrap_or("").to_string()),
                            label,
                            name,
                            product: attribute(attrs, "product").map(percent_decode),
                            accession: nucleotide_accession(attrs),
                        },
                    );
                }
            }
        }

        // Transcripts are the features that have exons or CDS blocks.
        let mut per_gene: HashMap<String, Vec<Transcript>> = HashMap::new();
        for id in exons.keys().chain(cds.keys()).collect::<std::collections::BTreeSet<_>>() {
            let Some(feature) = features.get(id) else { continue };
            let gene_feature = feature.parent.as_ref().and_then(|p| features.get(p)).filter(|g| g.kind != "mRNA");
            let gene = gene_feature.map(|g| g.label.clone()).unwrap_or_else(|| feature.label.clone());
            let gene_name = gene_feature.and_then(|g| g.name.clone()).or_else(|| feature.name.clone());
            let product = feature.product.clone().or_else(|| cds_products.get(id).cloned());
            let mut ex = exons.get(id).cloned().unwrap_or_default();
            let mut cd = cds.get(id).cloned().unwrap_or_default();
            if ex.is_empty() {
                ex = cd.clone(); // a transcript given only by its CDS
            }
            ex.sort();
            cd.sort();
            per_gene.entry(gene.clone()).or_default().push(Transcript {
                id: id.clone(),
                gene,
                contig: feature.contig.clone(),
                strand: feature.strand,
                exons: ex,
                cds: cd,
                transl_table: tables.get(id).copied().unwrap_or(1),
                gene_name,
                product,
                accession: feature.accession.clone(),
            });
        }
        let mut by_contig: HashMap<String, Vec<Transcript>> = HashMap::new();
        let mut transcripts_per_gene = HashMap::new();
        for (gene, mut txs) in per_gene {
            transcripts_per_gene.insert(gene, txs.len());
            // Canonical: the longest CDS, then the longest exon total, then the first id.
            txs.sort_by(|a, b| {
                let exon_len = |t: &Transcript| t.exons.iter().map(|e| e.1 - e.0 + 1).sum::<u64>();
                b.cds_len().cmp(&a.cds_len()).then(exon_len(b).cmp(&exon_len(a))).then(a.id.cmp(&b.id))
            });
            let canonical = txs.remove(0);
            by_contig.entry(canonical.contig.clone()).or_default().push(canonical);
        }
        for txs in by_contig.values_mut() {
            txs.sort_by(|a, b| a.start().cmp(&b.start()).then(a.id.cmp(&b.id)));
        }
        Ok(Annotation { by_contig, transcripts_per_gene })
    }
}

/// The reference sequence by contig name.
pub struct Genome(HashMap<String, Vec<u8>>);

impl Genome {
    pub fn new(contigs: Vec<(String, Vec<u8>)>) -> Genome {
        Genome(contigs.into_iter().collect())
    }

    pub(crate) fn base(&self, contig: &str, pos: u64) -> Option<u8> {
        self.0.get(contig)?.get((pos as usize).checked_sub(1)?).copied()
    }

    /// The length of a contig, or None when the reference has no such contig.
    pub(crate) fn len(&self, contig: &str) -> Option<u64> {
        self.0.get(contig).map(|s| s.len() as u64)
    }
}

/// What a variant does to the canonical transcript it falls in.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Effect {
    pub label: String,
    pub gene: Option<String>,
    pub transcript: Option<String>,
    pub aa: Option<String>,
    /// For a frameshift with a new stop: the length of the new protein and of the reference protein, as `new/reference` residues.
    pub truncation: Option<String>,
    pub n_transcripts: usize,
    pub gene_name: Option<String>,
    pub product: Option<String>,
    pub accession: Option<String>,
}

pub(crate) fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

/// Translate one codon. Table 1 is the standard code. Table 4 is the mold mitochondrial code
/// (TGA codes for W). Other tables fall back to the standard code.
fn translate(codon: [u8; 3], table: u8) -> u8 {
    const AA: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    let idx = |b: u8| match b {
        b'T' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    };
    match (idx(codon[0]), idx(codon[1]), idx(codon[2])) {
        (Some(a), Some(b), Some(c)) => {
            if table == 4 && codon == *b"TGA" {
                b'W'
            } else {
                AA[a * 16 + b * 4 + c]
            }
        }
        _ => b'X',
    }
}

fn effect_of(label: &str, tx: Option<&Transcript>, aa: Option<String>, ann: &Annotation) -> Effect {
    Effect {
        label: label.to_string(),
        gene: tx.map(|t| t.gene.clone()),
        transcript: tx.map(|t| t.id.clone()),
        aa,
        truncation: None,
        n_transcripts: tx.map(|t| ann.transcripts_per_gene.get(&t.gene).copied().unwrap_or(1)).unwrap_or(0),
        gene_name: tx.and_then(|t| t.gene_name.clone()),
        product: tx.and_then(|t| t.product.clone()),
        accession: tx.and_then(|t| t.accession.clone()),
    }
}

/// Where a position lies in a transcript.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Region {
    Cds,
    FivePrimeUtr,
    ThreePrimeUtr,
    NonCodingExon,
    Splice,
    Intron,
    Outside,
}

fn region_of(tx: &Transcript, pos: u64) -> Region {
    if tx.coding_index(pos).is_some() {
        return Region::Cds;
    }
    if tx.exons.iter().any(|&(s, e)| pos >= s && pos <= e) {
        // An exon base that is not in the CDS: a UTR, or a base of a non-coding transcript.
        let (Some(cds_start), Some(cds_end)) = (tx.cds.first().map(|c| c.0), tx.cds.last().map(|c| c.1)) else {
            return Region::NonCodingExon;
        };
        if pos >= cds_start && pos <= cds_end {
            return Region::NonCodingExon;
        }
        let low_side = pos < cds_start;
        let five_prime = (low_side && tx.strand == Strand::Plus) || (!low_side && tx.strand == Strand::Minus);
        return if five_prime { Region::FivePrimeUtr } else { Region::ThreePrimeUtr };
    }
    // Between exons: an intron derived from the gap. The first and last two bases are splice sites.
    for pair in tx.exons.windows(2) {
        let (intron_start, intron_end) = (pair[0].1 + 1, pair[1].0 - 1);
        if pos >= intron_start && pos <= intron_end {
            let splice = pos < intron_start + 2 || pos + 2 > intron_end;
            return if splice { Region::Splice } else { Region::Intron };
        }
    }
    Region::Outside
}

fn region_label(region: Region) -> &'static str {
    match region {
        Region::FivePrimeUtr => "5UTR_derived",
        Region::ThreePrimeUtr => "3UTR_derived",
        Region::NonCodingExon => "non_coding",
        Region::Splice => "splice_site",
        Region::Intron => "intron_derived",
        Region::Cds | Region::Outside => "unclassified",
    }
}

/// How serious an effect is, for choosing among transcripts that overlap a variant. Lower is more serious.
fn severity(label: &str) -> u8 {
    match label {
        "nonsense" => 0,
        "start_lost" => 1,
        "stop_lost" => 2,
        "frameshift" => 3,
        "splice_site" => 4,
        "large_indel_in_cds" => 5,
        "inframe_deletion" => 6,
        "inframe_insertion" => 7,
        "missense" => 8,
        "complex" => 9,
        "synonymous" => 10,
        "5UTR_derived" => 11,
        "3UTR_derived" => 12,
        "non_coding" => 13,
        "intron_derived" => 14,
        _ => 15,
    }
}

/// Classify one variant given in VCF form (an indel has one anchor base). Every canonical transcript
/// that overlaps the bases the variant affects is judged, and the most serious effect is returned.
pub fn classify_effect(ann: &Annotation, genome: &Genome, contig: &str, pos: u64, reference: &str, alt: &str) -> Effect {
    let Some(txs) = ann.by_contig.get(contig) else {
        return effect_of("intergenic", None, None, ann);
    };
    // The bases a variant affects: a SNP its base, a deletion the deleted bases (not the anchor), an
    // insertion the two bases on either side.
    let (span_start, span_end) = if reference.len() > alt.len() {
        let anchored = reference.as_bytes()[0] == alt.as_bytes()[0];
        let first = if anchored { pos + 1 } else { pos };
        (first, pos + reference.len() as u64 - 1)
    } else if reference.len() == alt.len() {
        (pos, pos + reference.len() as u64 - 1)
    } else {
        (pos, pos + 1)
    };
    let mut best: Option<Effect> = None;
    for tx in txs {
        if tx.start() > span_end || tx.end() < span_start {
            continue;
        }
        let effect = if reference.len() == 1 && alt.len() == 1 {
            snp_effect(ann, genome, tx, pos, alt.as_bytes()[0])
        } else if reference.len() == alt.len() {
            effect_of("unclassified", Some(tx), None, ann) // a multi-base substitution
        } else {
            indel_effect(ann, genome, tx, pos, reference, alt)
        };
        if best.as_ref().map_or(true, |b| severity(&effect.label) < severity(&b.label)) {
            best = Some(effect);
        }
    }
    best.unwrap_or_else(|| effect_of("intergenic", None, None, ann))
}

/// The effect of a SNP at `pos`, which lies between the first and last base of `tx`.
fn snp_effect(ann: &Annotation, genome: &Genome, tx: &Transcript, pos: u64, alt: u8) -> Effect {
    if let Some(idx) = tx.coding_index(pos) {
        return coding_snp(ann, genome, tx, idx, alt);
    }
    effect_of(region_label(region_of(tx, pos)), Some(tx), None, ann)
}

/// The effect of an insertion or deletion. The bases it affects are the deleted bases, or for an
/// insertion the two bases on either side. A splice site base anywhere makes it a splice site
/// variant. Otherwise all affected bases must lie in one region, or the variant is complex.
fn indel_effect(ann: &Annotation, genome: &Genome, tx: &Transcript, pos: u64, reference: &str, alt: &str) -> Effect {
    let symbolic = alt.starts_with('<');
    let net = reference.len() as i64 - alt.len() as i64;
    let affected: Vec<u64> = if !symbolic && net > 0 {
        let anchored = reference.as_bytes()[0] == alt.as_bytes()[0];
        let (first, count) = if anchored { (pos + 1, net as u64) } else { (pos, reference.len() as u64) };
        (first..first + count).collect()
    } else {
        vec![pos, pos + 1]
    };
    let regions: Vec<Region> = affected.iter().map(|&p| region_of(tx, p)).collect();
    let label = if regions.contains(&Region::Splice) {
        "splice_site"
    } else if regions.iter().all(|&r| r == regions[0]) {
        match regions[0] {
            Region::Cds if symbolic => "large_indel_in_cds",
            Region::Cds if net % 3 != 0 => {
                // The protein change and the truncation of a frameshift come from translating the shifted frame.
                let mut e = effect_of("frameshift", Some(tx), None, ann);
                if let Some((aa, truncation)) = frameshift_protein(genome, tx, pos, reference, alt) {
                    e.aa = Some(aa);
                    e.truncation = truncation;
                }
                return e;
            }
            Region::Cds => inframe_protein_effect(genome, tx, pos, reference, alt)
                .unwrap_or(if net > 0 { "inframe_deletion" } else { "inframe_insertion" }),
            other => region_label(other),
        }
    } else {
        "complex"
    };
    effect_of(label, Some(tx), None, ann)
}

fn coding_snp(ann: &Annotation, genome: &Genome, tx: &Transcript, idx: u64, alt: u8) -> Effect {
    let codon_no = idx / 3;
    let offset = (idx % 3) as usize;
    let mut ref_codon = [0u8; 3];
    for (k, slot) in ref_codon.iter_mut().enumerate() {
        let Some(pos) = tx.position_of(codon_no * 3 + k as u64) else {
            return effect_of("unclassified", Some(tx), None, ann);
        };
        let Some(base) = genome.base(&tx.contig, pos) else {
            return effect_of("unclassified", Some(tx), None, ann);
        };
        *slot = if tx.strand == Strand::Minus { complement(base.to_ascii_uppercase()) } else { base.to_ascii_uppercase() };
    }
    let mut alt_codon = ref_codon;
    alt_codon[offset] = if tx.strand == Strand::Minus { complement(alt.to_ascii_uppercase()) } else { alt.to_ascii_uppercase() };
    let (ref_aa, alt_aa) = (translate(ref_codon, tx.transl_table), translate(alt_codon, tx.transl_table));
    let label = if codon_no == 0 && ref_aa == b'M' && alt_aa != b'M' {
        "start_lost"
    } else if alt_aa == b'*' && ref_aa != b'*' {
        "nonsense"
    } else if ref_aa == b'*' && alt_aa != b'*' {
        "stop_lost"
    } else if ref_aa == alt_aa {
        "synonymous"
    } else {
        "missense"
    };
    let aa = format!("{}{}{}", ref_aa as char, codon_no + 1, alt_aa as char);
    effect_of(label, Some(tx), Some(aa), ann)
}

const INFO_DEFINITIONS: [&str; 7] = [
    "##INFO=<ID=EFFECT,Number=A,Type=String,Description=\"Effect on the canonical transcript, one per alternate allele (the _derived labels come from exons and CDS, not from UTR or intron features)\">",
    "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene (locus tag) of the canonical transcript\">",
    "##INFO=<ID=NUCCORE,Number=1,Type=String,Description=\"NCBI nucleotide accession of the canonical transcript, when the GFF3 gives one on the mRNA line (a GenBank entry in Dbxref, or a transcript_id that is an accession); used for a link to https://www.ncbi.nlm.nih.gov/nuccore/\">",
    "##INFO=<ID=TRANSCRIPT,Number=1,Type=String,Description=\"Canonical transcript: the one with the longest CDS\">",
    "##INFO=<ID=AA,Number=A,Type=String,Description=\"Protein change of a coding SNP, for example L424I; * is a stop. For a frameshift: the first changed residue, what it becomes, fs, and *N for the residues up to the new stop, for example E3Dfs*4 (E3Dfs*? when no stop is found)\">",
    "##INFO=<ID=TRUNC,Number=A,Type=String,Description=\"For a frameshift with a new stop: residues in the new protein / residues in the reference protein\">",
    "##INFO=<ID=NTX,Number=1,Type=Integer,Description=\"Number of transcripts of the gene\">",
];

/// Read a VCF and write it with `EFFECT`, `GENE`, `TRANSCRIPT`, `AA` and `NTX` added to INFO.
/// With no annotation every record is `EFFECT=unannotated`. The records keep their order and fields.
pub fn annotate_vcf<R: BufRead, W: io::Write>(
    reader: R,
    out: &mut W,
    ann: Option<&Annotation>,
    genome: &Genome,
) -> io::Result<()> {
    for line in reader.lines() {
        let line = line?;
        if line.starts_with("#CHROM") {
            for definition in INFO_DEFINITIONS {
                writeln!(out, "{definition}")?;
            }
            writeln!(out, "{line}")?;
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            writeln!(out, "{line}")?;
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, format!("VCF line has {} fields, needs at least 8: {line}", f.len())));
        }
        let pos: u64 = f[1]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("bad VCF position in: {line}")))?;
        let tags = match ann {
            None => "EFFECT=unannotated".to_string(),
            Some(ann) => {
                let effects: Vec<Effect> = f[4].split(',').map(|alt| classify_effect(ann, genome, f[0], pos, f[3], alt)).collect();
                let mut tags = format!("EFFECT={}", effects.iter().map(|e| e.label.as_str()).collect::<Vec<_>>().join(","));
                if let Some(first) = effects.first() {
                    if let (Some(gene), Some(tx)) = (&first.gene, &first.transcript) {
                        tags.push_str(&format!(";GENE={gene};TRANSCRIPT={tx}"));
                        if let Some(name) = &first.gene_name {
                            tags.push_str(&format!(";GENENAME={}", percent_encode_info(name)));
                        }
                        if let Some(product) = &first.product {
                            tags.push_str(&format!(";PRODUCT={}", percent_encode_info(product)));
                        }
                        if let Some(accession) = &first.accession {
                            tags.push_str(&format!(";NUCCORE={accession}"));
                        }
                    }
                    if effects.iter().any(|e| e.aa.is_some()) {
                        let aas: Vec<&str> = effects.iter().map(|e| e.aa.as_deref().unwrap_or(".")).collect();
                        tags.push_str(&format!(";AA={}", aas.join(",")));
                    }
                    if effects.iter().any(|e| e.truncation.is_some()) {
                        let t: Vec<&str> = effects.iter().map(|e| e.truncation.as_deref().unwrap_or(".")).collect();
                        tags.push_str(&format!(";TRUNC={}", t.join(",")));
                    }
                    if first.gene.is_some() {
                        tags.push_str(&format!(";NTX={}", first.n_transcripts));
                    }
                }
                tags
            }
        };
        let info = if f[7] == "." || f[7].is_empty() { tags } else { format!("{};{}", f[7], tags) };
        let mut fields: Vec<String> = f.iter().map(|s| s.to_string()).collect();
        fields[7] = info;
        writeln!(out, "{}", fields.join("\t"))?;
    }
    Ok(())
}

/// Annotate VCF text. `reference_fasta` is plain FASTA bytes and `gff` is plain GFF3 bytes. With no
/// GFF every record is `EFFECT=unannotated`. The browser build calls this; it does no file access.
pub fn annotate_text(vcf: &str, reference_fasta: &[u8], gff: Option<&[u8]>) -> Result<String, String> {
    let contigs = crate::read_fasta(reference_fasta).map_err(|e| format!("reference FASTA: {e}"))?;
    if contigs.is_empty() {
        return Err("the reference FASTA has no sequences".to_string());
    }
    let genome = Genome::new(contigs);
    let annotation = match gff {
        Some(bytes) => Some(Annotation::from_gff(bytes).map_err(|e| format!("GFF: {e}"))?),
        None => None,
    };
    let mut out = Vec::new();
    annotate_vcf(vcf.as_bytes(), &mut out, annotation.as_ref(), &genome).map_err(|e| format!("VCF: {e}"))?;
    String::from_utf8(out).map_err(|e| e.to_string())
}


fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

pub(crate) fn protein(codons: &[u8], table: u8) -> Vec<u8> {
    codons.chunks_exact(3).map(|c| translate([c[0], c[1], c[2]], table)).collect()
}

/// For an in-frame indel inside a CDS: does it change the start codon, create a premature stop, or
/// remove the terminal stop? Returns the label, or None when the protein changes in length only.
fn inframe_protein_effect(genome: &Genome, tx: &Transcript, pos: u64, reference: &str, alt: &str) -> Option<&'static str> {
    let cds = tx.coding_sequence(genome)?;
    let alt_cds = altered_cds(tx, &cds, pos, reference, alt)?;
    let (ref_aa, alt_aa) = (protein(&cds, tx.transl_table), protein(&alt_cds, tx.transl_table));
    if alt_aa.first() != ref_aa.first() && ref_aa.first() == Some(&b'M') {
        return Some("start_lost");
    }
    let ref_has_stop = ref_aa.last() == Some(&b'*');
    match alt_aa.iter().position(|&a| a == b'*') {
        Some(i) if i + 1 != alt_aa.len() || !ref_has_stop => Some("nonsense"),
        None if ref_has_stop => Some("stop_lost"),
        _ => None,
    }
}

/// The coding sequence with the indel applied. `reference` and `alt` are in VCF form (one anchor base).
pub(crate) fn altered_cds(tx: &Transcript, cds: &[u8], pos: u64, reference: &str, alt: &str) -> Option<Vec<u8>> {
    Some(if reference.len() > alt.len() {
        // The deleted bases follow the anchor. Their coding indices are contiguous.
        let deleted = (reference.len() - alt.len()) as u64;
        let indices: Vec<u64> = (pos + 1..=pos + deleted).map(|p| tx.coding_index(p)).collect::<Option<_>>()?;
        let (lo, hi) = (*indices.iter().min()? as usize, *indices.iter().max()? as usize + 1);
        [&cds[..lo], &cds[hi..]].concat()
    } else {
        let inserted = &alt.as_bytes()[reference.len()..];
        let anchor = tx.coding_index(pos)? as usize;
        let (at, bases) = match tx.strand {
            Strand::Plus => (anchor + 1, inserted.to_vec()),
            Strand::Minus => (anchor, reverse_complement(inserted)),
        };
        [&cds[..at], &bases[..], &cds[at..]].concat()
    })
}

/// The exon bases after the CDS, in transcript order (the 3' UTR, spliced). A frameshift near the end of the CDS can read into it.
pub(crate) fn three_prime_bases(genome: &Genome, tx: &Transcript) -> Vec<u8> {
    let (Some(first), Some(last)) = (tx.cds.first(), tx.cds.last()) else { return Vec::new() };
    let mut out = Vec::new();
    match tx.strand {
        Strand::Plus => {
            for &(s, e) in &tx.exons {
                for p in s.max(last.1 + 1)..=e {
                    if let Some(b) = genome.base(&tx.contig, p) {
                        out.push(b.to_ascii_uppercase());
                    }
                }
            }
        }
        Strand::Minus => {
            for &(s, e) in tx.exons.iter().rev() {
                if s >= first.0 {
                    continue;
                }
                for p in (s..=e.min(first.0 - 1)).rev() {
                    if let Some(b) = genome.base(&tx.contig, p) {
                        out.push(complement(b.to_ascii_uppercase()));
                    }
                }
            }
        }
    }
    out
}

/// What a frameshift does to the protein. Returns the protein change and, when the shifted frame has a stop, the truncation.
/// The change is the first residue that differs from the reference, what it becomes, `fs`, and `*N`: N residues from the first
/// changed one up to and including the new stop (`*?` when no stop is found, even in the 3' UTR). A first changed codon that is
/// already a stop is written `E3*fs`. The truncation is `new length/reference length` in residues.
fn frameshift_protein(genome: &Genome, tx: &Transcript, pos: u64, reference: &str, alt: &str) -> Option<(String, Option<String>)> {
    let cds = tx.coding_sequence(genome)?;
    let mut alt_cds = altered_cds(tx, &cds, pos, reference, alt)?;
    alt_cds.extend(three_prime_bases(genome, tx));
    let (ref_aa, alt_aa) = (protein(&cds, tx.transl_table), protein(&alt_cds, tx.transl_table));
    let ref_len = ref_aa.len() - usize::from(ref_aa.last() == Some(&b'*'));
    let common = ref_aa.len().min(alt_aa.len());
    let at = (0..common).find(|&i| ref_aa[i] != alt_aa[i]).unwrap_or(common);
    let ref_char = ref_aa.get(at).copied().unwrap_or(b'?') as char;
    let stop = alt_aa.iter().position(|&a| a == b'*').filter(|&j| j >= at);
    Some(match (alt_aa.get(at).copied(), stop) {
        (Some(b'*'), _) => (format!("{ref_char}{}*fs", at + 1), Some(format!("{at}/{ref_len}"))),
        (Some(a), Some(j)) => (format!("{ref_char}{}{}fs*{}", at + 1, a as char, j - at + 1), Some(format!("{j}/{ref_len}"))),
        (Some(a), None) => (format!("{ref_char}{}{}fs*?", at + 1, a as char), None),
        (None, _) => (format!("{ref_char}{}fs*?", at + 1), None),
    })
}

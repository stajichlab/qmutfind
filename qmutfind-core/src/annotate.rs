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
    fn coding_index(&self, pos: u64) -> Option<u64> {
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

struct Feature {
    kind: String,
    contig: String,
    strand: Strand,
    parent: Option<String>,
    label: String,
}

impl Annotation {
    /// Read a GFF3 file.
    pub fn from_gff<R: BufRead>(reader: R) -> io::Result<Annotation> {
        let mut features: HashMap<String, Feature> = HashMap::new();
        let mut exons: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let mut cds: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
        let mut tables: HashMap<String, u8> = HashMap::new();
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
                        }
                    }
                }
                kind => {
                    let Some(id) = attribute(attrs, "ID") else { continue };
                    let label = attribute(attrs, "locus_tag")
                        .or_else(|| attribute(attrs, "Name"))
                        .unwrap_or(id)
                        .to_string();
                    features.insert(
                        id.to_string(),
                        Feature {
                            kind: kind.to_string(),
                            contig: f[0].to_string(),
                            strand,
                            parent: attribute(attrs, "Parent").map(|p| p.split(',').next().unwrap_or("").to_string()),
                            label,
                        },
                    );
                }
            }
        }

        // Transcripts are the features that have exons or CDS blocks.
        let mut per_gene: HashMap<String, Vec<Transcript>> = HashMap::new();
        for id in exons.keys().chain(cds.keys()).collect::<std::collections::BTreeSet<_>>() {
            let Some(feature) = features.get(id) else { continue };
            let gene = feature
                .parent
                .as_ref()
                .and_then(|p| features.get(p))
                .filter(|g| g.kind != "mRNA")
                .map(|g| g.label.clone())
                .unwrap_or_else(|| feature.label.clone());
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

    fn base(&self, contig: &str, pos: u64) -> Option<u8> {
        self.0.get(contig)?.get((pos as usize).checked_sub(1)?).copied()
    }
}

/// What a variant does to the canonical transcript it falls in.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Effect {
    pub label: String,
    pub gene: Option<String>,
    pub transcript: Option<String>,
    pub aa: Option<String>,
    pub n_transcripts: usize,
}

fn complement(b: u8) -> u8 {
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
        n_transcripts: tx.map(|t| ann.transcripts_per_gene.get(&t.gene).copied().unwrap_or(1)).unwrap_or(0),
    }
}

/// Classify one variant given in VCF form (an indel has one anchor base).
pub fn classify_effect(ann: &Annotation, genome: &Genome, contig: &str, pos: u64, reference: &str, alt: &str) -> Effect {
    let none = effect_of("unclassified", None, None, ann);
    let Some(txs) = ann.by_contig.get(contig) else { return none };
    for tx in txs {
        if pos < tx.start() || pos > tx.end() {
            continue;
        }
        if reference.len() == 1 && alt.len() == 1 {
            if let Some(idx) = tx.coding_index(pos) {
                return coding_snp(ann, genome, tx, idx, alt.as_bytes()[0]);
            }
        }
    }
    none
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

//! The sequence around one variant, as aligned text rows: a position ruler, a marker under the changed columns, the reference DNA, its
//! protein, the variant DNA and its protein. Exon bases are capitals and all other bases small letters. The rows are on the gene's strand
//! (a minus-strand gene is shown reverse-complemented, so that the protein reads from left to right); outside a gene they are on the forward
//! strand. An insertion adds a column with a gap `-` in the reference row, a deletion puts `-` in the variant row.
//!
//! A protein letter sits under the middle base of its codon. Where a codon is split by an intron, the letter sits under the middle base of
//! the codon in the spliced coding sequence, and there are no letters in the intron.

use crate::annotate::{altered_cds, classify_effect, complement, protein, three_prime_bases, Strand};
use crate::{Annotation, Genome, Transcript};


/// The rows of the view. Every row has one character for each column.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SequenceContext {
    /// One line: the variant, the gene and the effect.
    pub header: String,
    /// (label, characters), in order: position, marker, reference, ref protein, variant, var protein. The protein rows are missing outside a gene.
    pub rows: Vec<(String, String)>,
    /// For each column: it is a changed column (a changed base, a deleted base or an inserted base).
    pub changed: Vec<bool>,
    /// For each column: the variant protein has a different letter from the reference protein.
    pub aa_changed: Vec<bool>,
    /// Remarks under the rows.
    pub notes: Vec<String>,
}

const LABEL_WIDTH: usize = 12;

struct Col {
    ref_base: Option<u8>,
    alt_base: Option<u8>,
    gpos: Option<u64>,
    changed: bool,
    in_exon: bool,
    /// Index of the reference base along the coding sequence (transcript direction), when it is in the CDS.
    ref_cds: Option<u64>,
    /// For an inserted column: it lies inside the CDS.
    alt_coding: bool,
}

impl SequenceContext {
    /// The view as text: the header, one line per row (label, then the characters), then the notes.
    pub fn to_text(&self) -> String {
        let mut out = format!("{}\n", self.header);
        for (label, chars) in &self.rows {
            out.push_str(&format!("{label:<LABEL_WIDTH$}{chars}\n"));
        }
        for note in &self.notes {
            out.push_str(&format!("note: {note}\n"));
        }
        out
    }

    /// The same text for a web page (inside `<pre>`): the changed columns of the marker and the DNA rows are in `<span class="ctx-var">`, the changed
    /// letters of the variant protein in `<span class="ctx-aa">`. Everything else is escaped text.
    pub fn to_html(&self) -> String {
        let esc = |c: char, out: &mut String| match c {
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '&' => out.push_str("&amp;"),
            other => out.push(other),
        };
        let mut out = String::new();
        self.header.chars().for_each(|c| esc(c, &mut out));
        out.push('\n');
        for (label, chars) in &self.rows {
            label.chars().for_each(|c| esc(c, &mut out));
            out.push_str(&" ".repeat(LABEL_WIDTH.saturating_sub(label.chars().count())));
            let marks: Option<(&[bool], &str)> = match label.as_str() {
                "marker" | "reference" | "variant" => Some((&self.changed, "ctx-var")),
                "var protein" => Some((&self.aa_changed, "ctx-aa")),
                _ => None,
            };
            let mut open = false;
            for (i, c) in chars.chars().enumerate() {
                let on = marks.map_or(false, |(flags, _)| flags.get(i).copied().unwrap_or(false));
                if on && !open {
                    out.push_str(&format!("<span class=\"{}\">", marks.map(|m| m.1).unwrap_or("")));
                    open = true;
                } else if !on && open {
                    out.push_str("</span>");
                    open = false;
                }
                esc(c, &mut out);
            }
            if open {
                out.push_str("</span>");
            }
            out.push('\n');
        }
        for note in &self.notes {
            out.push_str("note: ");
            note.chars().for_each(|c| esc(c, &mut out));
            out.push('\n');
        }
        out
    }
}

/// The view of the variant `contig:pos reference>alt` (VCF form: an indel has an anchor base) with `flank` bases of the reference on each side
/// (at least 1). Genes come from `ann`. Errors: the contig is not in the reference, the reference allele does not match it, or the variant is
/// not a substitution, an insertion or a deletion in VCF form.
pub fn sequence_context(genome: &Genome, ann: &Annotation, contig: &str, pos: u64, reference: &str, alt: &str, flank: u64) -> Result<SequenceContext, String> {
    let contig_len = genome.len(contig).ok_or_else(|| format!("the reference has no contig {contig}"))?;
    let (reference, alt) = (reference.to_ascii_uppercase(), alt.to_ascii_uppercase());
    if reference.is_empty() || alt.is_empty() || pos == 0 || pos + reference.len() as u64 - 1 > contig_len {
        return Err(format!("{contig}:{pos} is outside the reference"));
    }
    for (i, b) in reference.bytes().enumerate() {
        let have = genome.base(contig, pos + i as u64).map(|b| b.to_ascii_uppercase());
        if have != Some(b) {
            let seen: String = (0..reference.len() as u64).filter_map(|i| genome.base(contig, pos + i)).map(|b| b.to_ascii_uppercase() as char).collect();
            return Err(format!("the reference allele {reference} does not match the reference sequence at {contig}:{pos} ({seen})"));
        }
    }
    let flank = flank.max(1);
    let ref_len = reference.len() as u64;
    let kind = if reference.len() == alt.len() {
        Kind::Substitution
    } else if reference.len() > alt.len() && reference.starts_with(&alt) {
        Kind::Deletion
    } else if alt.len() > reference.len() && alt.starts_with(&reference) {
        Kind::Insertion
    } else {
        return Err("this kind of variant (a replacement of unequal length) is not drawn".to_string());
    };

    // The gene to show: the one that the annotation assigns the variant to, else a canonical transcript in the window.
    let effect = classify_effect(ann, genome, contig, pos, &reference, &alt);
    let (lo, hi) = (pos.saturating_sub(flank).max(1), (pos + ref_len - 1 + flank).min(contig_len));
    let in_variant = ann.transcripts_in(contig, pos, pos + ref_len - 1);
    let tx: Option<&Transcript> = effect
        .transcript
        .as_ref()
        .and_then(|id| in_variant.iter().find(|t| &t.id == id).copied())
        .or_else(|| ann.transcripts_in(contig, lo, hi).into_iter().next());

    // The columns, in genomic order.
    let exon_at = |p: u64| tx.map_or(false, |t| t.exons.iter().any(|&(s, e)| p >= s && p <= e));
    let cds_at = |p: u64| tx.and_then(|t| t.coding_index(p));
    let mut cols: Vec<Col> = Vec::new();
    for p in lo..=hi {
        let base = genome.base(contig, p).map(|b| b.to_ascii_uppercase()).unwrap_or(b'N');
        let off = p as i64 - pos as i64;
        let (alt_base, changed) = match kind {
            Kind::Substitution if off >= 0 && (off as usize) < alt.len() => (Some(alt.as_bytes()[off as usize]), alt.as_bytes()[off as usize] != base),
            Kind::Deletion if off >= alt.len() as i64 && off < ref_len as i64 => (None, true),
            _ => (Some(base), false),
        };
        cols.push(Col { ref_base: Some(base), alt_base, gpos: Some(p), changed, in_exon: exon_at(p), ref_cds: cds_at(p), alt_coding: false });
        if let Kind::Insertion = kind {
            if p == pos + ref_len - 1 {
                // The inserted bases follow the anchor; they are exon (coding) bases when the bases on both sides are.
                let both = |f: &dyn Fn(u64) -> bool| f(p) && (p + 1 > contig_len || f(p + 1));
                let (exon, coding) = (both(&|q| exon_at(q)), both(&|q| cds_at(q).is_some()));
                for b in alt.as_bytes()[reference.len()..].iter() {
                    cols.push(Col { ref_base: None, alt_base: Some(*b), gpos: None, changed: true, in_exon: exon, ref_cds: None, alt_coding: coding });
                }
            }
        }
    }
    let minus = tx.map_or(false, |t| t.strand == Strand::Minus);
    if minus {
        cols.reverse();
        for c in cols.iter_mut() {
            c.ref_base = c.ref_base.map(complement);
            c.alt_base = c.alt_base.map(complement);
        }
    }
    let n = cols.len();
    let letter = |b: Option<u8>, exon: bool| -> char {
        match b {
            None => '-',
            Some(b) if exon => b as char,
            Some(b) => (b as char).to_ascii_lowercase(),
        }
    };
    let ref_row: String = cols.iter().map(|c| letter(c.ref_base, c.in_exon)).collect();
    let alt_row: String = cols.iter().map(|c| letter(c.alt_base, c.in_exon)).collect();
    let marker: String = cols.iter().map(|c| if c.changed { '^' } else { ' ' }).collect();

    // The position ruler: the coordinate at every base whose position ends in 0, when it fits.
    let mut ruler = vec![' '; n];
    let mut next_free = 0usize;
    for (i, c) in cols.iter().enumerate() {
        let Some(p) = c.gpos else { continue };
        if p % 10 != 0 || i < next_free { continue }
        let label = p.to_string();
        if i + label.len() > n { continue }
        for (k, ch) in label.chars().enumerate() { ruler[i + k] = ch; }
        next_free = i + label.len() + 1;
    }

    let mut rows = vec![
        ("position".to_string(), ruler.into_iter().collect::<String>()),
        ("marker".to_string(), marker),
        ("reference".to_string(), ref_row.clone()),
    ];
    let mut notes = Vec::new();
    let mut aa_changed = vec![false; n];
    let coding_tx = tx.filter(|t| !t.cds.is_empty());
    if let Some(t) = coding_tx {
        if let Some(cds) = t.coding_sequence(genome) {
            let ref_prot: Vec<u8> = truncate_at_stop(protein(&cds, t.transl_table));
            // The CDS after the variant, or None when its protein cannot be given.
            let touches = |p: u64| t.coding_index(p).is_some();
            let span: Vec<u64> = match kind {
                Kind::Substitution => (pos..pos + ref_len).collect(),
                Kind::Deletion => (pos + alt.len() as u64..pos + ref_len).collect(),
                Kind::Insertion => vec![pos],
            };
            let any_coding = span.iter().any(|&p| touches(p));
            let all_coding = span.iter().all(|&p| touches(p));
            let alt_cds: Option<Vec<u8>> = if !any_coding {
                Some(cds.clone()) // in a UTR, an intron or outside the CDS: the protein is the same
            } else {
                match kind {
                    Kind::Substitution if all_coding => {
                        let mut a = cds.clone();
                        for (i, &p) in span.iter().enumerate() {
                            let idx = t.coding_index(p).unwrap_or(0) as usize;
                            let b = alt.as_bytes()[i];
                            a[idx] = if t.strand == Strand::Minus { complement(b) } else { b };
                        }
                        Some(a)
                    }
                    Kind::Deletion | Kind::Insertion if all_coding => altered_cds(t, &cds, pos, &reference, &alt),
                    _ => None,
                }
            };
            let ref_letters = protein_row(&cols, |c| c.ref_cds, &ref_prot);
            rows.push(("ref protein".to_string(), ref_letters.iter().collect()));
            rows.push(("variant".to_string(), alt_row.clone()));
            match alt_cds {
                Some(mut a) => {
                    if a.len() % 3 != cds.len() % 3 || a.len() != cds.len() {
                        a.extend(three_prime_bases(genome, t)); // a frameshift can read on into the 3' UTR
                    }
                    let alt_prot = truncate_at_stop(protein(&a, t.transl_table));
                    // The index of each variant-row base along the new coding sequence: the coding bases that are still there, in order.
                    let base_idx = cols.iter().find_map(|c| c.ref_cds);
                    let mut idxs: Vec<Option<u64>> = vec![None; n];
                    if let Some(first) = base_idx {
                        let mut count = 0u64;
                        for (i, c) in cols.iter().enumerate() {
                            let coding = c.ref_cds.is_some() || c.alt_coding;
                            if coding && c.alt_base.is_some() {
                                idxs[i] = Some(first + count);
                                count += 1;
                            }
                        }
                    }
                    let alt_letters = protein_row(&cols, |_| None, &[]).iter().enumerate().map(|(i, _)| {
                        match idxs[i] {
                            Some(k) if k % 3 == 1 => alt_prot.get((k / 3) as usize).map(|&b| b as char).unwrap_or(' '),
                            _ => ' ',
                        }
                    }).collect::<Vec<char>>();
                    for i in 0..n {
                        aa_changed[i] = alt_letters[i] != ' ' && alt_letters[i] != ref_letters[i];
                    }
                    rows.push(("var protein".to_string(), alt_letters.into_iter().collect()));
                }
                None => {
                    rows.push(("var protein".to_string(), " ".repeat(n)));
                    notes.push("the variant reaches the edge of the coding sequence, so its protein is not shown".to_string());
                }
            }
        } else {
            rows.push(("variant".to_string(), alt_row.clone()));
            notes.push("the coding sequence of this gene could not be read".to_string());
        }
    } else {
        rows.push(("variant".to_string(), alt_row.clone()));
        notes.push(if tx.is_some() { "the transcript has no coding sequence (non-coding RNA): no protein is shown".to_string() } else { "the variant is not in a gene: the rows show the forward strand and no protein".to_string() });
    }
    if tx.is_some() && effect.transcript.is_none() {
        notes.push("the variant is not in this gene; the gene is shown because it lies nearby".to_string());
    }
    notes.push("capitals are exon bases and small letters are not; ^ marks the changed columns; - is a gap".to_string());

    let aa = effect.aa.as_deref().filter(|a| *a != ".").map(|a| format!(" {a}")).unwrap_or_default();
    let header = match tx {
        Some(t) => format!(
            "{contig}:{pos} {reference}>{alt}  gene {} ({}, {} strand{})  effect {}{aa}",
            t.gene, t.id, if minus { "minus" } else { "plus" },
            if minus { "; the rows show the gene strand, the reverse complement of the reference" } else { "" },
            effect.label
        ),
        None => format!("{contig}:{pos} {reference}>{alt}  effect {}", effect.label),
    };
    let changed: Vec<bool> = cols.iter().map(|c| c.changed).collect();
    Ok(SequenceContext { header, rows, changed, aa_changed, notes })
}

enum Kind {
    Substitution,
    Deletion,
    Insertion,
}

/// The protein up to and including its first stop.
fn truncate_at_stop(mut prot: Vec<u8>) -> Vec<u8> {
    if let Some(i) = prot.iter().position(|&b| b == b'*') {
        prot.truncate(i + 1);
    }
    prot
}

/// One letter per column: the amino acid under the middle base of its codon (`index` gives the coding index of a column), a space elsewhere.
fn protein_row(cols: &[Col], index: impl Fn(&Col) -> Option<u64>, prot: &[u8]) -> Vec<char> {
    cols.iter()
        .map(|c| match index(c) {
            Some(i) if i % 3 == 1 => prot.get((i / 3) as usize).map(|&b| b as char).unwrap_or(' '),
            _ => ' ',
        })
        .collect()
}

/// The view as text or HTML from a reference FASTA and a GFF3 file (both plain bytes), like `annotate_text`.
pub fn context_text(reference_fasta: &[u8], gff: Option<&[u8]>, contig: &str, pos: u64, reference: &str, alt: &str, flank: u64, html: bool) -> Result<String, String> {
    let contigs = crate::read_fasta(reference_fasta).map_err(|e| format!("reference FASTA: {e}"))?;
    if contigs.is_empty() {
        return Err("the reference FASTA has no sequences".to_string());
    }
    let genome = Genome::new(contigs);
    let ann = match gff {
        Some(bytes) => Annotation::from_gff(bytes).map_err(|e| format!("GFF: {e}"))?,
        None => Annotation::empty(),
    };
    let view = sequence_context(&genome, &ann, contig, pos, reference, alt, flank)?;
    Ok(if html { view.to_html() } else { view.to_text() })
}

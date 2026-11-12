//! Left-alignment of indels. An indel in a repeat can be written at several positions.
//! The leftmost form is the standard one, and it is what `bcftools norm` writes.

use crate::Call;

/// Return the leftmost form of one indel, or `None` when it cannot be anchored to the left.
fn left_align(seq: &[u8], mut pos: u64, reference: &str, alt: &str) -> Option<(u64, String, String)> {
    let mut r: Vec<u8> = reference.bytes().collect();
    let mut a: Vec<u8> = alt.bytes().collect();
    loop {
        // Trim the bases the two alleles share at the right end.
        while !r.is_empty() && !a.is_empty() && r.last() == a.last() {
            r.pop();
            a.pop();
        }
        if !r.is_empty() && !a.is_empty() {
            break;
        }
        // One allele is empty: take the base to the left into both alleles.
        if pos <= 1 {
            return None;
        }
        pos -= 1;
        let base = *seq.get((pos - 1) as usize)?;
        r.insert(0, base);
        a.insert(0, base);
    }
    // Keep one anchor base: trim shared bases at the left while both alleles stay longer than one.
    while r.len() > 1 && a.len() > 1 && r[0] == a[0] {
        r.remove(0);
        a.remove(0);
        pos += 1;
    }
    Some((pos, String::from_utf8(r).ok()?, String::from_utf8(a).ok()?))
}

/// Write every explicit indel in its leftmost form. SNPs, symbolic calls and indels that
/// cannot be anchored to the left are returned as they were.
pub fn left_align_calls(calls: Vec<Call>, reference: &[(String, Vec<u8>)]) -> Vec<Call> {
    calls
        .into_iter()
        .map(|mut c| {
            if c.svlen.is_none() && c.reference.len() != c.alt.len() {
                if let Some((_, seq)) = reference.iter().find(|(name, _)| *name == c.contig) {
                    if let Some((pos, r, a)) = left_align(seq, c.pos, &c.reference, &c.alt) {
                        c.pos = pos;
                        c.reference = r;
                        c.alt = a;
                    }
                }
            }
            c
        })
        .collect()
}

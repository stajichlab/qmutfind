//! SAM input. A record is reduced to what the caller needs: where it starts, its CIGAR
//! operations without clips, and the aligned bases in reference orientation.

use crate::CigarOp;

/// One mapped SAM record. `ref_start` is 0-based. `seq` holds only the aligned bases:
/// soft-clipped bases are removed. SAM stores SEQ in reference orientation already.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SamRecord {
    pub contig: String,
    pub ref_start: usize,
    pub cigar_ops: Vec<CigarOp>,
    pub seq: Vec<u8>,
}

/// Parse one SAM line. Returns `Ok(None)` for header lines and for unmapped (0x4) and
/// secondary (0x100) records. Supplementary records (0x800) are kept. A malformed line is
/// an error, never a silent skip.
pub fn parse_sam_line(line: &str) -> Result<Option<SamRecord>, String> {
    if line.starts_with('@') {
        return Ok(None);
    }
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 11 {
        return Err(format!("SAM line has {} fields, needs at least 11: {line}", f.len()));
    }
    let flag: u32 = f[1].parse().map_err(|_| format!("bad SAM flag '{}'", f[1]))?;
    if flag & 0x4 != 0 || flag & 0x100 != 0 {
        return Ok(None);
    }
    let pos1: usize = f[3].parse().map_err(|_| format!("bad SAM position '{}'", f[3]))?;
    if pos1 == 0 {
        return Err(format!("mapped SAM record has position 0: {line}"));
    }

    let mut ops = Vec::new();
    let (mut lead_clip, mut trail_clip) = (0usize, 0usize);
    let mut digits = String::new();
    let mut seen_aligned_op = false;
    for ch in f[5].chars() {
        if ch.is_ascii_digit() {
            digits.push(ch);
            continue;
        }
        let n: u32 = digits.parse().map_err(|_| format!("bad CIGAR '{}'", f[5]))?;
        digits.clear();
        match ch {
            'M' | '=' | 'X' => {
                ops.push(CigarOp::Match(n));
                seen_aligned_op = true;
            }
            'I' => {
                ops.push(CigarOp::Insertion(n));
                seen_aligned_op = true;
            }
            'D' => {
                ops.push(CigarOp::Deletion(n));
                seen_aligned_op = true;
            }
            'S' if !seen_aligned_op => lead_clip += n as usize,
            'S' => trail_clip += n as usize,
            'H' => {}
            other => return Err(format!("unsupported CIGAR operation '{other}' in '{}'", f[5])),
        }
    }

    let full = f[9].as_bytes();
    if lead_clip + trail_clip > full.len() {
        return Err(format!("CIGAR clips exceed the sequence length: {line}"));
    }
    let seq: Vec<u8> = full[lead_clip..full.len() - trail_clip].iter().map(|b| b.to_ascii_uppercase()).collect();
    let consumed: usize = ops
        .iter()
        .map(|op| match op {
            CigarOp::Match(n) | CigarOp::Insertion(n) => *n as usize,
            CigarOp::Deletion(_) => 0,
        })
        .sum();
    if consumed != seq.len() {
        return Err(format!("CIGAR '{}' does not match the sequence length: {line}", f[5]));
    }
    Ok(Some(SamRecord { contig: f[2].to_string(), ref_start: pos1 - 1, cigar_ops: ops, seq }))
}

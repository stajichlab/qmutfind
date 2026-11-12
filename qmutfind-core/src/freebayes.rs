//! Read a FreeBayes VCF into calls. A multi-base substitution is split into SNPs, as
//! `vcfallelicprimitives` does, so this step needs no external program. Records are kept
//! with the same depth and fraction thresholds as the built-in caller.

use std::io::{self, BufRead};

use crate::{Call, CallerConfig};

fn invalid(message: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

fn info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='))
}

fn only_acgt(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

/// Read calls from a FreeBayes VCF. A record needs `DP` and `AO` in INFO.
pub fn calls_from_freebayes_vcf<R: BufRead>(reader: R, config: &CallerConfig) -> io::Result<Vec<Call>> {
    let mut calls = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            return Err(invalid(format!("VCF line has {} fields, needs at least 8: {line}", f.len())));
        }
        let pos: u64 = f[1].parse().map_err(|_| invalid(format!("bad VCF position '{}'", f[1])))?;
        let depth: u32 = info_value(f[7], "DP")
            .ok_or_else(|| invalid(format!("record has no DP: {line}")))?
            .parse()
            .map_err(|_| invalid(format!("bad DP in: {line}")))?;
        let aos: Vec<u32> = info_value(f[7], "AO")
            .ok_or_else(|| invalid(format!("record has no AO: {line}")))?
            .split(',')
            .map(|v| v.parse().map_err(|_| invalid(format!("bad AO in: {line}"))))
            .collect::<io::Result<_>>()?;
        let reference = f[3].to_ascii_uppercase();
        for (alt, &alt_count) in f[4].split(',').zip(aos.iter()) {
            let alt = alt.to_ascii_uppercase();
            if depth < config.min_depth
                || depth == 0
                || (alt_count as f64 / depth as f64)
                    < if reference.len() == alt.len() { config.min_snp_fraction } else { config.min_indel_fraction }
                || !only_acgt(&reference)
                || !only_acgt(&alt)
            {
                continue;
            }
            let make = |p: u64, r: &str, a: &str| Call {
                contig: f[0].to_string(),
                pos: p,
                reference: r.to_string(),
                alt: a.to_string(),
                depth,
                alt_count,
                svlen: None,
            };
            if reference.len() == alt.len() && reference.len() > 1 {
                // A multi-base substitution: one SNP for each differing base.
                for (i, (r, a)) in reference.bytes().zip(alt.bytes()).enumerate() {
                    if r != a {
                        calls.push(make(pos + i as u64, &(r as char).to_string(), &(a as char).to_string()));
                    }
                }
            } else {
                calls.push(make(pos, &reference, &alt));
            }
        }
    }
    Ok(calls)
}

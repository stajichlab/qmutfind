//! VCF output. One record per sample per mutation, with the sample name in INFO.

use std::io::{self, Write};

use crate::{classify, Call, VariantClass};

fn class_name(class: VariantClass) -> &'static str {
    match class {
        VariantClass::Snp => "SNP",
        VariantClass::SmallIndel => "SMALL_INDEL",
        VariantClass::LargeIndel => "LARGE_INDEL",
    }
}

/// Write `calls` as VCF 4.2. `contigs` holds (name, length) in reference order.
/// `large_indel_min` is the size cutoff for the large indel class.
pub fn write_vcf<W: Write>(
    out: &mut W,
    contigs: &[(String, usize)],
    sample: &str,
    calls: &[Call],
    large_indel_min: usize,
) -> io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    for (name, len) in contigs {
        writeln!(out, "##contig=<ID={name},length={len}>")?;
    }
    writeln!(out, "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample label\">")?;
    writeln!(out, "##INFO=<ID=CLASS,Number=1,Type=String,Description=\"SNP, SMALL_INDEL or LARGE_INDEL\">")?;
    writeln!(out, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth at the site (at the anchor base for an indel)\">")?;
    writeln!(out, "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Reads that support the alternate allele\">")?;
    writeln!(out, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Alternate allele fraction, AC divided by DP\">")?;
    writeln!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
    for c in calls {
        let class = class_name(classify(&c.reference, &c.alt, large_indel_min));
        let af = c.alt_count as f64 / c.depth as f64;
        writeln!(
            out,
            "{}\t{}\t.\t{}\t{}\t.\tPASS\tSAMPLE={};CLASS={};DP={};AC={};AF={:.4}",
            c.contig, c.pos, c.reference, c.alt, sample, class, c.depth, c.alt_count, af
        )?;
    }
    Ok(())
}

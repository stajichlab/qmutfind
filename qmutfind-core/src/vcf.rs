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
    writeln!(out, "##INFO=<ID=AO,Number=A,Type=Integer,Description=\"Reads that support the alternate allele\">")?;
    writeln!(out, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele fraction, AO divided by DP\">")?;
    writeln!(out, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">")?;
    writeln!(out, "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Length of the structural variant, estimated from reads\">")?;
    writeln!(out, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Position and length are estimates from reads\">")?;
    writeln!(out, "##ALT=<ID=INS,Description=\"Insertion whose sequence is not given\">")?;
    writeln!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
    for c in calls {
        let class = match c.svlen {
            Some(len) if len >= large_indel_min => VariantClass::LargeIndel,
            Some(_) => VariantClass::SmallIndel,
            None => classify(&c.reference, &c.alt, large_indel_min),
        };
        let af = c.alt_count as f64 / c.depth as f64;
        let sv = match c.svlen {
            Some(len) => format!(";SVTYPE=INS;SVLEN={len};IMPRECISE"),
            None => String::new(),
        };
        writeln!(
            out,
            "{}\t{}\t.\t{}\t{}\t.\tPASS\tSAMPLE={};CLASS={};DP={};AO={};AF={:.4}{}",
            c.contig, c.pos, c.reference, c.alt, sample, class_name(class), c.depth, c.alt_count, af, sv
        )?;
    }
    Ok(())
}

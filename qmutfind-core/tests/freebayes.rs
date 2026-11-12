//! Reading a FreeBayes VCF into calls, with the same depth and fraction thresholds as the built-in caller.

use qmutfind_core::{calls_from_freebayes_vcf, CallerConfig};

fn config() -> CallerConfig {
    CallerConfig { min_depth: 5, min_alt_fraction: 0.5, large_indel_min: 50 }
}

/// One FreeBayes-style record. `ao` is the AO value as written (comma-separated for several alleles).
fn record(pos: u64, reference: &str, alt: &str, dp: u32, ao: &str) -> String {
    format!("chrA\t{pos}\t.\t{reference}\t{alt}\t50\t.\tAB=0;AO={ao};DP={dp};RO=1\tGT:DP\t1:{dp}\n")
}

fn parse(text: &str) -> Vec<qmutfind_core::Call> {
    let vcf = format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\n{text}");
    calls_from_freebayes_vcf(vcf.as_bytes(), &config()).unwrap()
}

#[test]
fn a_snp_record_becomes_one_call_with_depth_and_alt_count() {
    let calls = parse(&record(100, "A", "G", 20, "18"));

    assert_eq!(calls.len(), 1);
    let c = &calls[0];
    assert_eq!((c.contig.as_str(), c.pos, c.reference.as_str(), c.alt.as_str()), ("chrA", 100, "A", "G"));
    assert_eq!((c.depth, c.alt_count, c.svlen), (20, 18, None));
}

#[test]
fn records_below_the_depth_or_fraction_thresholds_are_dropped() {
    let calls = parse(&format!(
        "{}{}{}",
        record(100, "A", "G", 20, "18"), // kept
        record(200, "A", "G", 20, "4"),  // fraction 0.2
        record(300, "A", "G", 4, "4"),   // depth 4 is below 5
    ));

    assert_eq!(calls.iter().map(|c| c.pos).collect::<Vec<_>>(), vec![100]);
}

#[test]
fn a_multi_base_substitution_is_split_into_snps_at_the_differing_positions_only() {
    let calls = parse(&record(100, "ACG", "TCA", 20, "19"));

    let snps: Vec<_> = calls.iter().map(|c| (c.pos, c.reference.as_str(), c.alt.as_str())).collect();
    assert_eq!(snps, vec![(100, "A", "T"), (102, "G", "A")]);
    assert!(calls.iter().all(|c| c.depth == 20 && c.alt_count == 19));
}

#[test]
fn each_allele_of_a_multi_allele_record_uses_its_own_ao_value() {
    let calls = parse(&record(100, "A", "G,T", 12, "8,3"));

    assert_eq!(calls.len(), 1); // G has 8 of 12, T has 3 of 12
    assert_eq!((calls[0].alt.as_str(), calls[0].alt_count), ("G", 8));
}

#[test]
fn an_indel_record_passes_through_unchanged() {
    let calls = parse(&record(100, "CAA", "C", 20, "17"));

    assert_eq!(calls.len(), 1);
    assert_eq!((calls[0].reference.as_str(), calls[0].alt.as_str()), ("CAA", "C"));
}

#[test]
fn a_record_without_dp_is_an_error() {
    let vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrA\t100\t.\tA\tG\t50\t.\tAO=9\n";

    assert!(calls_from_freebayes_vcf(vcf.as_bytes(), &config()).is_err());
}

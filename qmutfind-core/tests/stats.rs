//! Summary statistics of a variant set, from a VCF. The fixture is small and checked by hand.

use qmutfind_core::compute_stats;

const VCF: &str = "##fileformat=VCFv4.2
##contig=<ID=c1,length=100000>
##contig=<ID=c2,length=50000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
c1\t100\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1.0000;EFFECT=missense
c1\t150\t.\tC\tT\t.\tPASS\tCLASS=SNP;AF=0.9000;EFFECT=synonymous
c1\t5000\t.\tA\tC\t.\tPASS\tCLASS=SNP;AF=0.5000;EFFECT=intergenic
c1\t12000\t.\tG\tT\t.\tPASS\tCLASS=SNP;AF=0.3000;EFFECT=nonsense
c1\t12050\t.\tCAA\tC\t.\tPASS\tCLASS=SMALL_INDEL;AF=0.7000;EFFECT=frameshift
c2\t10\t.\tT\tTGGG\t.\tPASS\tCLASS=SMALL_INDEL;AF=1.0000;EFFECT=inframe_insertion
c2\t20000\t.\tA\tG\t.\tPASS\tAF=0.8500
";

fn stats() -> qmutfind_core::Stats {
    compute_stats(VCF.as_bytes(), 50).unwrap()
}

#[test]
fn counts_are_given_by_class_and_by_effect_and_a_missing_class_is_derived_from_the_alleles() {
    let s = stats();

    assert_eq!(s.total, 7);
    assert_eq!(s.by_class["SNP"], 5); // the last record has no CLASS and is a SNP by its alleles
    assert_eq!(s.by_class["SMALL_INDEL"], 2);
    assert_eq!(s.by_class.get("LARGE_INDEL").copied().unwrap_or(0), 0);
    assert_eq!(s.by_effect["missense"], 1);
    assert_eq!(s.by_effect["frameshift"], 1);
    assert_eq!(s.by_effect["unannotated"], 1); // the record without an EFFECT tag
}

#[test]
fn ts_tv_has_its_counts_and_is_undefined_when_there_are_no_transversions() {
    let s = stats();

    assert_eq!((s.transitions, s.transversions), (3, 2));
    assert!((s.ts_tv.unwrap() - 1.5).abs() < 1e-9);

    let only_ts = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nc1\t5\t.\tA\tG\t.\t.\t.\n";
    assert_eq!(compute_stats(only_ts.as_bytes(), 50).unwrap().ts_tv, None);
}

#[test]
fn the_substitution_spectrum_has_twelve_types_and_six_strand_collapsed_classes() {
    let s = stats();

    assert_eq!(s.spectrum["A>G"], 2);
    assert_eq!(s.spectrum["C>T"], 1);
    assert_eq!(s.spectrum["A>C"], 1);
    assert_eq!(s.spectrum["G>T"], 1);
    // Collapsed to the pyrimidine of the pair: A>G is T>C, A>C is T>G, G>T is C>A.
    assert_eq!(s.spectrum_collapsed["T>C"], 2);
    assert_eq!(s.spectrum_collapsed["C>T"], 1);
    assert_eq!(s.spectrum_collapsed["T>G"], 1);
    assert_eq!(s.spectrum_collapsed["C>A"], 1);
}

#[test]
fn indel_lengths_are_signed_and_the_frameshift_fraction_counts_lengths_not_a_multiple_of_three() {
    let s = stats();

    assert_eq!(s.indel_lengths[&-2], 1); // a 2 bp deletion
    assert_eq!(s.indel_lengths[&3], 1); // a 3 bp insertion
    assert!((s.frameshift_length_fraction.unwrap() - 0.5).abs() < 1e-9);
}

#[test]
fn counts_and_density_per_contig_use_the_contig_lengths_from_the_header() {
    let s = stats();

    assert_eq!(s.per_contig["c1"].count, 5);
    assert_eq!(s.per_contig["c2"].count, 2);
    assert!((s.per_contig["c1"].per_kb - 0.05).abs() < 1e-9);
    assert!((s.per_contig["c2"].per_kb - 0.04).abs() < 1e-9);
    assert!((s.overall_per_kb.unwrap() - 7.0 / 150.0).abs() < 1e-9);
}

#[test]
fn the_fraction_histogram_has_ten_bins_with_one_in_the_last_and_mid_range_calls_are_flagged() {
    let s = stats();

    assert_eq!(s.fraction_histogram, vec![0, 0, 0, 1, 0, 1, 0, 1, 1, 3]);
    // Fractions from 0.2 to 0.8: 0.5, 0.3 and 0.7 of the 7 records with a fraction.
    assert_eq!(s.mid_fraction_calls, 3);
    assert!((s.mid_fraction_share.unwrap() - 3.0 / 7.0).abs() < 1e-9);
}

#[test]
fn clustering_counts_the_variants_with_a_neighbour_within_the_window_and_the_smallest_gap() {
    let s = stats(); // the window argument is 50 bp in these tests

    assert_eq!(s.smallest_gap, Some(50));
    // c1: 100 and 150 are 50 apart, and 12000 and 12050 are 50 apart.
    assert_eq!(s.clustered, 4);
}

#[test]
fn a_record_with_too_few_fields_is_an_error() {
    assert!(compute_stats("c1\t5\t.\tA\n".as_bytes(), 50).is_err());
}


#[test]
fn the_json_is_well_formed_and_holds_the_values_with_null_for_an_undefined_ratio() {
    let json: serde_json::Value = serde_json::from_str(&stats().to_json()).unwrap();

    assert_eq!(json["total"], 7);
    assert_eq!(json["by_class"]["SNP"], 5);
    assert_eq!(json["by_effect"]["unannotated"], 1);
    assert!((json["ts_tv"].as_f64().unwrap() - 1.5).abs() < 1e-9);
    assert_eq!(json["spectrum"]["A>G"], 2);
    assert_eq!(json["spectrum_collapsed"]["T>C"], 2);
    assert_eq!(json["indel_lengths"]["-2"], 1);
    assert_eq!(json["indel_lengths"]["3"], 1);
    assert_eq!(json["per_contig"]["c1"]["count"], 5);
    assert_eq!(json["fraction_histogram"], serde_json::json!([0, 0, 0, 1, 0, 1, 0, 1, 1, 3]));
    assert_eq!(json["mid_fraction_calls"], 3);
    assert_eq!(json["smallest_gap"], 50);
    assert_eq!(json["clustered"], 4);

    let no_tv = compute_stats("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nc1\t5\t.\tA\tG\t.\t.\t.\n".as_bytes(), 50).unwrap();
    let json: serde_json::Value = serde_json::from_str(&no_tv.to_json()).unwrap();
    assert!(json["ts_tv"].is_null());
}

#[test]
fn the_tsv_has_metric_key_value_columns_with_a_header() {
    let tsv = stats().to_tsv();
    let lines: Vec<&str> = tsv.lines().collect();

    assert_eq!(lines[0], "metric\tkey\tvalue");
    assert!(lines.contains(&"total\t\t7"));
    assert!(lines.contains(&"ts_tv\t\t1.5000"));
    assert!(lines.contains(&"transitions\t\t3"));
    assert!(lines.contains(&"spectrum\tA>G\t2"));
    assert!(lines.contains(&"by_effect\tmissense\t1"));
    assert!(lines.contains(&"indel_length\t-2\t1"));
    assert!(lines.contains(&"per_contig_per_kb\tc1\t0.0500"));
    assert!(lines.contains(&"fraction_histogram\t0.9-1.0\t3"));
    assert!(lines.contains(&"clustered\t\t4"));
}

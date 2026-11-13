//! Windowed read depth (from the caller's counts), its mosdepth-style BED text, and the genome-wide depth plot.

mod common;

use common::{assert_label_rules, contig_labels, crowded_contigs, scaffold_contigs};
use qmutfind_core::{depth_plot_svg, Caller, CallerConfig, CigarOp, SamRecord};

const LEN: usize = 1000;

fn reference() -> Vec<(String, Vec<u8>)> {
    vec![("c".to_string(), (0..LEN).map(|i| b"ACGT"[(i * 7 + i / 3) % 4]).collect())]
}

fn config() -> CallerConfig {
    CallerConfig { min_depth: 1, min_snp_fraction: 0.5, min_indel_fraction: 0.5, large_indel_min: 50 }
}

/// A read that copies the reference over `[start, start + len)`.
fn read(refs: &[(String, Vec<u8>)], start: usize, len: usize) -> SamRecord {
    SamRecord { contig: "c".to_string(), ref_start: start, cigar_ops: vec![CigarOp::Match(len as u32)], seq: refs[0].1[start..start + len].to_vec(), mapq: 60 }
}

#[test]
fn window_means_include_uncovered_bases_and_the_last_window_may_be_short() {
    let refs = reference();
    let mut caller = Caller::new(refs.clone(), config());
    for _ in 0..5 {
        caller.add_sam(&read(&refs, 0, 500)); // depth 5 on the first 500 bases only
    }

    let windows = caller.depth_windows(300); // 0-300, 300-600, 600-900, 900-1000

    let means: Vec<f64> = windows.iter().map(|w| w.mean).collect();
    assert_eq!(windows.len(), 4);
    assert!((means[0] - 5.0).abs() < 1e-9);
    assert!((means[1] - 5.0 * 200.0 / 300.0).abs() < 1e-9, "{means:?}"); // 200 of 300 bases covered
    assert_eq!(means[2], 0.0);
    assert_eq!((windows[3].start, windows[3].end), (900, 1000));
}

#[test]
fn the_bed_text_has_four_columns_like_mosdepth_regions() {
    let refs = reference();
    let mut caller = Caller::new(refs.clone(), config());
    caller.add_sam(&read(&refs, 0, 1000));

    let bed = caller.depth_bed(500);

    assert_eq!(bed, "c\t0\t500\t1.00\nc\t500\t1000\t1.00\n");
}

fn bed(rows: &[(&str, u64, u64, f64)]) -> String {
    rows.iter().map(|(c, s, e, d)| format!("{c}\t{s}\t{e}\t{d}\n")).collect()
}

/// The value of attribute `name` in the `n`-th element that contains `marker`.
fn attr(svg: &str, marker: &str, n: usize, name: &str) -> String {
    let start = svg.match_indices(marker).nth(n).unwrap_or_else(|| panic!("no element {n} with {marker}")).0;
    let element = &svg[start..start + svg[start..].find('>').unwrap()];
    let key = format!("{name}=\"");
    let from = element.find(&key).unwrap_or_else(|| panic!("no {name} in {element}")) + key.len();
    element[from..from + element[from..].find('"').unwrap()].to_string()
}

fn num(svg: &str, marker: &str, n: usize, name: &str) -> f64 {
    attr(svg, marker, n, name).parse().unwrap()
}

#[test]
fn a_deeper_window_is_drawn_higher() {
    let svg = depth_plot_svg(&bed(&[("a", 0, 1000, 10.0), ("a", 1000, 2000, 20.0)]), None, "input", None).unwrap();

    let low = num(&svg, "data-depth=\"10.00\"", 0, "data-y");
    let high = num(&svg, "data-depth=\"20.00\"", 0, "data-y");
    assert!(high < low, "{high} {low}");
}

#[test]
fn contigs_can_be_ordered_by_input_by_length_or_by_natural_name() {
    let rows = bed(&[("chr10", 0, 500, 5.0), ("chr2", 0, 3000, 5.0), ("chr1", 0, 1000, 5.0)]);
    let order = |how: &str| -> Vec<String> {
        let svg = depth_plot_svg(&rows, None, how, None).unwrap();
        svg.match_indices("data-contig-label=\"").map(|(i, _)| {
            let rest = &svg[i + "data-contig-label=\"".len()..];
            rest[..rest.find('"').unwrap()].to_string()
        }).collect()
    };

    assert_eq!(order("input"), vec!["chr10", "chr2", "chr1"]);
    assert_eq!(order("length"), vec!["chr2", "chr1", "chr10"]);
    assert_eq!(order("name"), vec!["chr1", "chr2", "chr10"]);
}

#[test]
fn an_unknown_order_is_an_error_that_lists_the_choices() {
    let err = depth_plot_svg(&bed(&[("a", 0, 100, 5.0)]), None, "sideways", None).unwrap_err();

    assert!(err.contains("input") && err.contains("length") && err.contains("name"), "{err}");
}

#[test]
fn a_window_far_above_the_median_is_clipped_and_marked() {
    let rows = bed(&[("a", 0, 1000, 10.0), ("a", 1000, 2000, 10.0), ("a", 2000, 3000, 10.0), ("a", 3000, 4000, 500.0)]);

    let svg = depth_plot_svg(&rows, None, "input", None).unwrap();

    assert!(svg.contains("data-clipped=\"true\""));
    assert!(svg.contains("data-median=\"10.00\""));
    assert_eq!(svg.matches("data-clipped=\"true\"").count(), 1);
}

#[test]
fn rows_of_a_run_length_bed_are_averaged_by_their_length_into_each_bin() {
    // One contig of 1,000 bp as a per-base style BED: 100 bp at depth 100, then 900 bp at depth 0, drawn in one bin.
    let svg = depth_plot_svg(&bed(&[("a", 0, 100, 100.0), ("a", 100, 1000, 0.0)]), Some(&[("a".to_string(), 1000)]), "input", None);
    let svg = svg.unwrap();

    let depth: f64 = attr(&svg, "class=\"depth\"", 0, "data-depth").parse().unwrap();
    assert!((depth - 10.0).abs() < 1e-6, "{depth}");
}

#[test]
fn contig_lengths_default_to_the_largest_end_in_the_bed_and_a_given_list_overrides_them() {
    let rows = bed(&[("a", 0, 1000, 5.0)]);

    let svg = depth_plot_svg(&rows, Some(&[("a".to_string(), 1000), ("b".to_string(), 2000)]), "input", None).unwrap();

    assert!(svg.contains("data-contig-label=\"b\""));
}

#[test]
fn variant_ticks_are_drawn_when_a_vcf_is_given() {
    let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\na\t10\t.\tA\tG\t.\tPASS\tCLASS=SNP\na\t900\t.\tAC\tA\t.\tPASS\tCLASS=SMALL_INDEL\n";

    let svg = depth_plot_svg(&bed(&[("a", 0, 1000, 5.0)]), None, "input", Some(vcf)).unwrap();

    assert_eq!(svg.matches("class=\"vtick\"").count(), 2);
}

#[test]
fn an_empty_bed_is_an_error_and_a_bad_row_names_the_line() {
    assert!(depth_plot_svg("", None, "input", None).unwrap_err().contains("no rows"));
    let err = depth_plot_svg("a\t0\tx\t5\n", None, "input", None).unwrap_err();
    assert!(err.contains("line 1"), "{err}");
}

#[test]
fn a_bed_row_on_a_contig_missing_from_a_given_list_is_an_error() {
    let err = depth_plot_svg(&bed(&[("z", 0, 100, 5.0)]), Some(&[("a".to_string(), 100)]), "input", None).unwrap_err();

    assert!(err.contains("z"), "{err}");
}


// Contig labels: bold, turned 90 degrees, right-justified below the variant tick row, and never on top of each other.

fn depth_plot_of(contigs: &[(String, usize)], vcf: Option<&str>) -> String {
    let rows: String = contigs.iter().map(|(n, l)| format!("{n}\t0\t{}\t10.0\n", (*l).min(1000))).collect();
    depth_plot_svg(&rows, Some(contigs), "input", vcf).unwrap()
}

fn one_call_vcf(contig: &str) -> String {
    format!("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{contig}\t500\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1.0\n")
}

#[test]
fn the_labels_of_the_depth_plot_follow_the_label_rules_and_hang_below_the_variant_tick_row() {
    // The variant ticks end at y 262, so no label may start above that.
    let contigs = scaffold_contigs();
    assert_label_rules(&depth_plot_of(&contigs, Some(&one_call_vcf(&contigs[0].0))), 262.0);
    assert_label_rules(&depth_plot_of(&contigs, None), 262.0);
}

#[test]
fn in_the_depth_plot_too_the_wider_contig_keeps_its_label_when_two_collide() {
    let labels = contig_labels(&depth_plot_of(&crowded_contigs(), None));

    let shown = |name: &str| !labels.iter().find(|l| l.name == name).unwrap().text.is_empty();
    assert!(shown("wider") && !shown("narrow"));
}

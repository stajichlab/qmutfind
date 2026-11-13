//! The whole-genome plot: where the variants are along the contigs, drawn as SVG.

mod common;

use common::{assert_label_rules, contig_labels, crowded_contigs, scaffold_contigs, viewbox_height};
use qmutfind_core::genome_plot_svg;

const HEADER: &str = "##fileformat=VCFv4.2
##contig=<ID=chrA,length=1000>
##contig=<ID=chrB,length=3000>
##contig=<ID=chrC,length=500>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
";

fn vcf(records: &[&str]) -> String {
    let mut text = HEADER.to_string();
    for r in records {
        text.push_str(r);
        text.push('\n');
    }
    text
}

const CALLS: [&str; 4] = [
    "chrA\t100\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1.0000",
    "chrA\t400\t.\tCT\tC\t.\tPASS\tCLASS=SMALL_INDEL;AF=0.5000",
    "chrB\t2000\t.\tG\t<INS>\t.\tPASS\tCLASS=LARGE_INDEL;AF=0.9000;SVLEN=120",
    "chrB\t2500\t.\tT\tC\t.\tPASS\tCLASS=SNP;AF=0.8000",
];

/// The value of attribute `name` in the `n`-th element that contains `marker`.
fn attr_of(svg: &str, marker: &str, n: usize, name: &str) -> String {
    let start = svg.match_indices(marker).nth(n).unwrap_or_else(|| panic!("no element {n} with {marker}")).0;
    let end = start + svg[start..].find('>').unwrap();
    let element = &svg[start..end];
    let key = format!("{name}=\"");
    let from = element.find(&key).unwrap_or_else(|| panic!("no {name} in {element}")) + key.len();
    element[from..from + element[from..].find('"').unwrap()].to_string()
}

fn num(svg: &str, marker: &str, n: usize, name: &str) -> f64 {
    attr_of(svg, marker, n, name).parse().unwrap()
}

#[test]
fn there_is_one_marker_per_call_with_the_right_shape_for_its_class() {
    let svg = genome_plot_svg(&vcf(&CALLS), None, "input").unwrap();

    assert_eq!(svg.matches("class=\"call\"").count(), 4);
    assert_eq!(svg.matches("data-class=\"SNP\"").count(), 2);
    assert_eq!(svg.matches("data-class=\"SMALL_INDEL\"").count(), 1);
    assert_eq!(svg.matches("data-class=\"LARGE_INDEL\"").count(), 1);
    assert!(svg.contains("<circle class=\"call\" data-class=\"SNP\""));
    assert!(svg.contains("<rect class=\"call\" data-class=\"SMALL_INDEL\""));
    assert!(svg.contains("<polygon class=\"call\" data-class=\"LARGE_INDEL\""));
}

#[test]
fn every_contig_is_labelled_even_when_it_has_no_calls_and_calls_follow_the_contig_order() {
    let svg = genome_plot_svg(&vcf(&CALLS), None, "input").unwrap();

    for name in ["chrA", "chrB", "chrC"] {
        assert!(svg.contains(&format!("data-contig-label=\"{name}\"")), "{name} is not labelled");
    }
    let a = num(&svg, "data-contig=\"chrA\"", 0, "data-x");
    let b = num(&svg, "data-contig=\"chrB\"", 0, "data-x");
    assert!(b > a);
}

#[test]
fn a_position_maps_to_x_in_proportion_to_the_genome_length() {
    // chrA is 1000 bp and chrB 3000 bp, so chrB is three times as wide. Calls at the same fraction of each contig
    // (10% in) are 100 bp into chrA and 300 bp into chrB.
    let svg = genome_plot_svg(
        &vcf(&["chrA\t101\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1", "chrB\t301\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1"]),
        None,
        "input",
    )
    .unwrap();

    let band_a = num(&svg, "data-band=\"chrA\"", 0, "width");
    let band_b = num(&svg, "data-band=\"chrB\"", 0, "width");
    assert!((band_b / band_a - 3.0).abs() < 0.01, "{band_a} {band_b}");
    let start_a = num(&svg, "data-band=\"chrA\"", 0, "x");
    let start_b = num(&svg, "data-band=\"chrB\"", 0, "x");
    let offset_a = num(&svg, "data-contig=\"chrA\"", 0, "data-x") - start_a;
    let offset_b = num(&svg, "data-contig=\"chrB\"", 0, "data-x") - start_b;
    assert!((offset_a / band_a - 0.1).abs() < 0.002, "{offset_a} {band_a}");
    assert!((offset_b / band_b - 0.1).abs() < 0.002, "{offset_b} {band_b}");
}

#[test]
fn a_higher_allele_fraction_is_drawn_higher_on_the_page() {
    let svg = genome_plot_svg(&vcf(&CALLS), None, "input").unwrap();

    let y_af_1 = num(&svg, "data-pos=\"100\"", 0, "data-y");
    let y_af_half = num(&svg, "data-pos=\"400\"", 0, "data-y");
    assert!(y_af_1 < y_af_half, "{y_af_1} {y_af_half}");
}

#[test]
fn the_density_bins_add_up_to_the_number_of_calls() {
    let svg = genome_plot_svg(&vcf(&CALLS), None, "input").unwrap();

    let total: usize = svg
        .match_indices("class=\"bin\"")
        .map(|(i, _)| attr_of(&svg[i..], "class=\"bin\"", 0, "data-count").parse::<usize>().unwrap())
        .sum();
    assert_eq!(total, 4);
}

#[test]
fn contig_lengths_can_come_from_the_reference_when_the_vcf_has_no_contig_lines() {
    let bare = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrA\t100\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1\n";
    let contigs = vec![("chrA".to_string(), 1000usize), ("chrB".to_string(), 2000usize)];

    let svg = genome_plot_svg(bare, Some(&contigs), "input").unwrap();

    assert!(svg.contains("data-contig-label=\"chrB\""));
}

#[test]
fn without_any_contig_lengths_the_error_says_what_to_give() {
    let bare = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrA\t100\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1\n";

    let err = genome_plot_svg(bare, None, "input").unwrap_err();

    assert!(err.contains("contig") && err.contains("reference"), "{err}");
}

#[test]
fn a_call_on_a_contig_with_no_known_length_is_an_error_not_a_silent_drop() {
    let err = genome_plot_svg(&vcf(&["chrZ\t5\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1"]), None, "input").unwrap_err();

    assert!(err.contains("chrZ"), "{err}");
}

#[test]
fn a_vcf_without_calls_still_gives_a_plot_that_says_there_are_none() {
    let svg = genome_plot_svg(&vcf(&[]), None, "input").unwrap();

    assert_eq!(svg.matches("class=\"call\"").count(), 0);
    assert!(svg.contains("No variants"));
}

#[test]
fn the_class_is_worked_out_from_the_alleles_when_the_class_tag_is_missing() {
    let svg = genome_plot_svg(
        &vcf(&["chrA\t10\t.\tA\tG\t.\tPASS\tAF=1", "chrA\t20\t.\tAC\tA\t.\tPASS\tAF=1", "chrA\t30\t.\tA\tAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\t.\tPASS\tAF=1"]),
        None,
        "input",
    )
    .unwrap();

    assert_eq!(svg.matches("data-class=\"SNP\"").count(), 1);
    assert_eq!(svg.matches("data-class=\"SMALL_INDEL\"").count(), 1);
    assert_eq!(svg.matches("data-class=\"LARGE_INDEL\"").count(), 1);
}

#[test]
fn a_call_without_an_allele_fraction_is_drawn_on_a_separate_lane_and_the_plot_says_so() {
    let svg = genome_plot_svg(&vcf(&["chrA\t10\t.\tA\tG\t.\tPASS\tCLASS=SNP"]), None, "input").unwrap();

    assert!(svg.contains("no AF in the VCF"));
}

#[test]
fn contig_names_are_escaped_and_the_output_is_deterministic() {
    let text = "##fileformat=VCFv4.2\n##contig=<ID=a<b&c,length=100>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\na<b&c\t5\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1\n";

    let first = genome_plot_svg(text, None, "input").unwrap();

    assert!(!first.contains("a<b&c"));
    assert!(first.contains("a&lt;b&amp;c"));
    assert_eq!(first, genome_plot_svg(text, None, "input").unwrap());
}

#[test]
fn calls_far_apart_on_one_contig_fall_in_different_bins_and_close_calls_share_one() {
    let svg = genome_plot_svg(
        &vcf(&[
            "chrB\t30\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1",
            "chrB\t31\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1",
            "chrB\t2970\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1",
        ]),
        None,
        "input",
    )
    .unwrap();

    let counts: Vec<usize> = svg
        .match_indices("class=\"bin\"")
        .map(|(i, _)| attr_of(&svg[i..], "class=\"bin\"", 0, "data-count").parse().unwrap())
        .collect();
    assert_eq!(counts, vec![2, 1], "bins in order along the contig");
    let first_x = num(&svg, "class=\"bin\"", 0, "x");
    let last_x = num(&svg, "class=\"bin\"", 1, "x");
    assert!(last_x - first_x > 500.0, "the two bins are at the two ends of chrB: {first_x} {last_x}");
}

#[test]
fn the_variant_plot_orders_contigs_by_input_length_or_natural_name_like_the_depth_plot() {
    let text = "##fileformat=VCFv4.2\n##contig=<ID=chr10,length=500>\n##contig=<ID=chr2,length=3000>\n##contig=<ID=chr1,length=1000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let order = |how: &str| -> Vec<String> {
        let svg = genome_plot_svg(text, None, how).unwrap();
        svg.match_indices("data-contig-label=\"").map(|(i, _)| {
            let rest = &svg[i + "data-contig-label=\"".len()..];
            rest[..rest.find('"').unwrap()].to_string()
        }).collect()
    };

    assert_eq!(order("input"), vec!["chr10", "chr2", "chr1"]);
    assert_eq!(order("length"), vec!["chr2", "chr1", "chr10"]);
    assert_eq!(order("name"), vec!["chr1", "chr2", "chr10"]);
    assert!(genome_plot_svg(text, None, "sideways").is_err());
}


// Contig labels: bold, turned 90 degrees, right-justified below the plot, and never on top of each other.

const EMPTY_VCF: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn plot_of(contigs: &[(String, usize)]) -> String {
    genome_plot_svg(EMPTY_VCF, Some(contigs), "input").unwrap()
}

#[test]
fn the_labels_of_the_variant_plot_follow_the_label_rules_and_hang_below_the_allele_fraction_area() {
    // The allele fraction area ends at y 300, so no label may start above it.
    assert_label_rules(&genome_plot_svg(&vcf(&CALLS), None, "input").unwrap(), 300.0);
    assert_label_rules(&plot_of(&scaffold_contigs()), 300.0);
}

#[test]
fn when_labels_would_collide_the_wider_contig_keeps_its_label() {
    let labels = contig_labels(&plot_of(&crowded_contigs()));

    let shown = |name: &str| !labels.iter().find(|l| l.name == name).unwrap().text.is_empty();
    assert!(shown("wider") && !shown("narrow"), "the narrow contig comes first in the file but the wider one gets the label");
}

#[test]
fn a_contig_without_a_label_still_names_itself_in_a_tooltip_on_its_band() {
    let svg = plot_of(&crowded_contigs());

    assert!(svg.contains("data-band=\"narrow\""));
    let band = &svg[svg.find("data-band=\"narrow\"").unwrap()..];
    assert!(band[..band.find("</rect>").unwrap()].contains("<title>narrow: 2300 bp</title>"));
}

#[test]
fn a_long_contig_name_is_cut_with_an_ellipsis_and_the_tooltip_keeps_the_full_name() {
    let name = "scaffold_with_a_very_long_name_that_does_not_fit_1234567890";
    let svg = plot_of(&[(name.to_string(), 5_000_000)]);

    let label = &contig_labels(&svg)[0];
    assert!(label.text.ends_with('…') && label.text.chars().count() <= 24, "{}", label.text);
    assert!(svg.contains(&format!("<title>{name}: 5000000 bp</title>")));
}

#[test]
fn the_legend_is_below_the_labels_and_inside_the_picture() {
    let svg = plot_of(&scaffold_contigs());

    let label_bottom = contig_labels(&svg).iter().filter(|l| !l.text.is_empty()).map(|l| l.y + l.text.chars().count() as f64 * 6.6).fold(0.0, f64::max);
    let legend_y = num(&svg, "<circle cx=\"74.0\"", 0, "cy");
    assert!(legend_y > label_bottom, "legend at {legend_y}, labels end at {label_bottom}");
    assert!(legend_y + 10.0 <= viewbox_height(&svg));
}

//! The locus view: gene models and mutations drawn as SVG. The gene model is the hand-checked
//! fixture of tests/effects.rs: contig c1, plus strand, exons 21-50 and 81-110, CDS 24-50 and 81-104.

use qmutfind_core::{parse_region, read_report_variants, render_locus_svg, render_report_html, Annotation};

const GFF: &str = "##gff-version 3
c1\tsrc\tgene\t21\t110\t.\t+\t.\tID=gene1;locus_tag=G1
c1\tsrc\tmRNA\t21\t110\t.\t+\t.\tID=tx1;Parent=gene1
c1\tsrc\texon\t21\t50\t.\t+\t.\tParent=tx1
c1\tsrc\texon\t81\t110\t.\t+\t.\tParent=tx1
c1\tsrc\tCDS\t24\t50\t.\t+\t0\tParent=tx1
c1\tsrc\tCDS\t81\t104\t.\t+\t0\tParent=tx1
c2\tsrc\tgene\t11\t22\t.\t-\t.\tID=gene2;locus_tag=G2
c2\tsrc\tmRNA\t11\t22\t.\t-\t.\tID=tx2;Parent=gene2
c2\tsrc\texon\t11\t22\t.\t-\t.\tParent=tx2
c2\tsrc\tCDS\t11\t22\t.\t-\t0\tParent=tx2
";

const VCF: &str = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
c1\t28\t.\tC\tA\t.\tPASS\tSAMPLE=s1;EFFECT=missense;GENE=G1;TRANSCRIPT=tx1;AA=A2D
c1\t60\t.\tA\tC\t.\tPASS\tSAMPLE=s1;EFFECT=intron_derived;GENE=G1;TRANSCRIPT=tx1
c1\t200\t.\tT\tA\t.\tPASS\tSAMPLE=s1;EFFECT=intergenic
c1\t999\t.\tT\tA\t.\tPASS\tSAMPLE=s1;EFFECT=missense
c1\t70\t.\tA\tC\t.\tPASS\tSAMPLE=s1
";

fn svg(region: &str) -> String {
    let ann = Annotation::from_gff(GFF.as_bytes()).unwrap();
    let variants = read_report_variants(VCF.as_bytes()).unwrap();
    render_locus_svg(&ann, &variants, &parse_region(region).unwrap())
}

fn count(haystack: &str, needle: &str) -> usize {
    haystack.matches(needle).count()
}

/// The value of attribute `name` in the first element that contains `marker`.
fn attr(svg: &str, marker: &str, name: &str) -> String {
    let element_start = svg.find(marker).unwrap_or_else(|| panic!("no element with {marker}"));
    let tail = &svg[element_start..];
    let key = format!("{name}=\"");
    let at = tail.find(&key).unwrap() + key.len();
    tail[at..tail[at..].find('"').map(|i| at + i).unwrap()].to_string()
}

#[test]
fn a_region_is_parsed_with_or_without_thousands_separators_and_bad_ones_are_rejected() {
    let r = parse_region("c1:21-110").unwrap();
    assert_eq!((r.contig.as_str(), r.start, r.end), ("c1", 21, 110));
    let r = parse_region("NC_007198.1:3,851,274-3,852,983").unwrap();
    assert_eq!((r.contig.as_str(), r.start, r.end), ("NC_007198.1", 3_851_274, 3_852_983));
    for bad in ["c1", "c1:10", "c1:20-10", "c1:a-b", "c1:0-5", ":1-5"] {
        assert!(parse_region(bad).is_err(), "{bad} should be rejected");
    }
}

#[test]
fn the_gene_model_has_cds_blocks_derived_utrs_an_intron_line_a_label_and_the_strand() {
    let s = svg("c1:1-200");

    assert!(s.starts_with("<svg") && s.contains("viewBox="));
    assert_eq!(count(&s, "class=\"cds\""), 2);
    assert_eq!(count(&s, "class=\"utr\""), 2); // 21-23 and 105-110, derived from the exons
    assert_eq!(count(&s, "class=\"intron\""), 1);
    assert!(s.contains(">G1<"));
    assert!(s.contains("data-strand=\"+\""));
    assert!(svg("c2:1-40").contains("data-strand=\"-\""));
}

#[test]
fn each_variant_inside_the_region_is_one_circle_with_its_effect_and_a_tooltip() {
    let s = svg("c1:1-200");

    assert_eq!(count(&s, "class=\"variant\""), 4); // 28, 60, 70 and 200 are in the region; 999 is not
    assert!(s.contains("data-effect=\"missense\""));
    assert!(s.contains("data-effect=\"intron_derived\""));
    assert!(s.contains("c1:28 C&gt;A missense A2D G1"), "{s}");
}

#[test]
fn a_variant_without_an_effect_tag_is_drawn_as_unannotated_and_one_outside_the_region_is_not_drawn() {
    let s = svg("c1:1-200");

    // Position 70 has no EFFECT tag. It is drawn with the unannotated colour. Position 999 is outside.
    assert!(s.contains("data-effect=\"unannotated\""));
    assert!(!s.contains("c1:999"));
    assert_eq!(count(&s, "class=\"variant\""), 4);
}

#[test]
fn the_x_position_follows_the_mapping_from_the_region_to_the_plot_width() {
    // The plot is 880 wide from x=60. A base is drawn at its centre: x = 60 + (pos - start + 0.5) / length * 880.
    let s = svg("c1:21-110");
    let cx: f64 = attr(&s, "data-pos=\"28\"", "cx").parse().unwrap();

    let expected = 60.0 + (28.0 - 21.0 + 0.5) / 90.0 * 880.0;
    assert!((cx - expected).abs() < 0.01, "cx {cx} expected {expected}");
}

#[test]
fn the_html_page_wraps_the_svg_supports_dark_mode_and_escapes_the_title() {
    let html = render_report_html("A <b>&\"title\"", "<svg></svg>", &[("calls".to_string(), "4".to_string())]);

    assert!(html.starts_with("<!doctype html>"));
    assert!(html.contains("<svg></svg>"));
    assert!(html.contains("prefers-color-scheme: dark"));
    assert!(html.contains("A &lt;b&gt;&amp;&quot;title&quot;"));
    assert!(!html.contains("<b>"));
    assert!(html.contains("calls") && html.contains(">4<"));
}

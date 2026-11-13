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

// ---- Impact labels next to each variant, and a table of all variants on the page ----

fn svg_for(vcf: &str, region: &str) -> String {
    let ann = Annotation::from_gff(GFF.as_bytes()).unwrap();
    let variants = read_report_variants(vcf.as_bytes()).unwrap();
    render_locus_svg(&ann, &variants, &parse_region(region).unwrap())
}

/// The text of the label element that has `data-pos="pos"`.
fn label_text(svg: &str, pos: u64) -> Option<String> {
    let marker = format!("class=\"variant-label\" data-pos=\"{pos}\"");
    let start = svg.find(&marker)?;
    let open_end = start + svg[start..].find('>')? + 1;
    Some(svg[open_end..open_end + svg[open_end..].find("</text>")?].to_string())
}

#[test]
fn a_variant_with_a_protein_change_is_labelled_with_it_next_to_its_circle() {
    let s = svg("c1:1-200");

    assert_eq!(label_text(&s, 28).as_deref(), Some("A2D"));
    let circle_x: f64 = attr(&s, "class=\"variant\" data-effect=\"missense\"", "cx").parse().unwrap();
    let label_x: f64 = attr(&s, "class=\"variant-label\" data-pos=\"28\"", "x").parse().unwrap();
    assert!(label_x > circle_x, "{label_x} {circle_x}");
    assert_eq!(attr(&s, "class=\"variant-label\" data-pos=\"28\"", "text-anchor"), "start");
}

#[test]
fn a_variant_without_a_protein_change_is_labelled_with_a_short_effect_name() {
    let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
c1\t60\t.\tA\tC\t.\tPASS\tEFFECT=intron_derived\n\
c1\t90\t.\tA\tC\t.\tPASS\tEFFECT=frameshift;GENE=G1\n\
c1\t150\t.\tA\tC\t.\tPASS\tEFFECT=5UTR_derived\n\
c1\t180\t.\tA\tC\t.\tPASS\t.\n";

    let s = svg_for(vcf, "c1:1-200");

    assert_eq!(label_text(&s, 60).as_deref(), Some("intron"));
    assert_eq!(label_text(&s, 90).as_deref(), Some("frameshift"));
    assert_eq!(label_text(&s, 150).as_deref(), Some("5UTR"));
    assert_eq!(label_text(&s, 180).as_deref(), Some("unannotated"));
}

#[test]
fn a_label_that_would_overlap_the_label_before_it_on_the_same_level_is_left_out_but_its_circle_stays() {
    // In 1,000 bases each base is 0.88 px wide, so positions 28 and 44 are 14 px apart: one level, but the labels are wider.
    let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
c1\t28\t.\tC\tA\t.\tPASS\tEFFECT=missense;AA=A2D\n\
c1\t44\t.\tC\tA\t.\tPASS\tEFFECT=missense;AA=L10P\n";

    let s = svg_for(vcf, "c1:1-1000");

    assert_eq!(count(&s, "class=\"variant\""), 2);
    // The first label would cross the second circle, which is on the same level, so only the second is drawn.
    assert_eq!(count(&s, "class=\"variant-label\""), 1);
    assert_eq!(label_text(&s, 44).as_deref(), Some("L10P"));
    assert_eq!(label_text(&s, 28), None);
}

#[test]
fn a_label_near_the_right_edge_is_drawn_to_the_left_of_its_circle() {
    let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
c1\t199\t.\tC\tA\t.\tPASS\tEFFECT=missense;AA=A2DELTAXXXX\n";

    let s = svg_for(vcf, "c1:1-200");

    assert_eq!(attr(&s, "class=\"variant-label\" data-pos=\"199\"", "text-anchor"), "end");
    let circle_x: f64 = attr(&s, "class=\"variant\"", "cx").parse().unwrap();
    let label_x: f64 = attr(&s, "class=\"variant-label\" data-pos=\"199\"", "x").parse().unwrap();
    assert!(label_x < circle_x);
}

#[test]
fn the_page_lists_every_variant_in_the_region_in_a_table_even_when_its_label_was_left_out() {
    let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
c1\t28\t.\tC\tA\t.\tPASS\tEFFECT=missense;GENE=G1;AA=A2D\n\
c1\t44\t.\tC\tA\t.\tPASS\tEFFECT=missense;GENE=G1;AA=L10P\n\
c1\t5000\t.\tC\tA\t.\tPASS\tEFFECT=missense;GENE=G1;AA=Q9R\n";

    let page = qmutfind_core::report_page(vcf, Some(GFF.as_bytes()), "c1:1-1000", "t").unwrap();

    assert!(page.contains("<table class=\"variants\">"));
    for needle in ["c1:28", "A2D", "c1:44", "L10P", "G1", "missense"] {
        assert!(page.contains(needle), "{needle}");
    }
    assert!(!page.contains("Q9R"), "a variant outside the region is not listed");
}

#[test]
fn the_first_and_last_axis_labels_are_anchored_inside_the_plot_so_they_are_not_cut_off() {
    let s = svg("c1:1-200");

    let anchors: Vec<String> = s
        .match_indices("font-size=\"11\" text-anchor=\"")
        .map(|(i, m)| {
            let rest = &s[i + m.len()..];
            rest[..rest.find('"').unwrap()].to_string()
        })
        .collect();

    assert_eq!(anchors, vec!["start", "middle", "end"]);
}

// Direction of transcription: an arrow under each gene model and chevrons along its line.

const LEFT: f64 = 60.0;
const PLOT: f64 = 880.0;

/// The points of the first element that contains `marker`, as (x, y) pairs.
fn points(svg: &str, marker: &str) -> Vec<(f64, f64)> {
    attr(svg, marker, "points")
        .split_whitespace()
        .map(|p| {
            let (x, y) = p.split_once(',').unwrap();
            (x.parse().unwrap(), y.parse().unwrap())
        })
        .collect()
}

fn tip_and_base(svg: &str) -> ((f64, f64), f64) {
    let p = points(svg, "class=\"arrowhead\"");
    assert_eq!(p.len(), 3, "an arrowhead is a triangle");
    // The tip is the point on the centre line; the other two share an x (the base).
    let tip = *p.iter().find(|q| p.iter().filter(|r| (r.0 - q.0).abs() < 0.05).count() == 1).unwrap();
    let base = p.iter().find(|q| (q.0 - tip.0).abs() > 0.05).unwrap().0;
    (tip, base)
}

#[test]
fn a_plus_strand_gene_has_an_arrowhead_at_its_right_end_that_points_right() {
    let s = svg("c1:1-200"); // gene 21-110
    let x_end = LEFT + 110.0 / 200.0 * PLOT; // the right edge of base 110

    let (tip, base) = tip_and_base(&s);

    assert!((tip.0 - x_end).abs() < 0.2, "tip {} vs gene end {x_end}", tip.0);
    assert!(base < tip.0, "the base is left of the tip, so the arrow points right");
    assert_eq!(attr(&s, "class=\"arrowhead\"", "data-strand"), "+");
}

#[test]
fn a_minus_strand_gene_has_an_arrowhead_at_its_left_end_that_points_left() {
    let s = svg("c2:1-40"); // gene 11-22 on the minus strand
    let x_start = LEFT + 10.0 / 40.0 * PLOT; // the left edge of base 11

    let (tip, base) = tip_and_base(&s);

    assert!((tip.0 - x_start).abs() < 0.2, "tip {} vs gene start {x_start}", tip.0);
    assert!(base > tip.0, "the base is right of the tip, so the arrow points left");
    assert_eq!(attr(&s, "class=\"arrowhead\"", "data-strand"), "-");
}

#[test]
fn when_the_end_of_transcription_is_outside_the_region_the_arrow_points_out_of_the_plot_edge() {
    let (tip, base) = tip_and_base(&svg("c1:1-60")); // the plus gene continues to 110, past the right edge
    assert!((tip.0 - (LEFT + PLOT)).abs() < 0.2 && base < tip.0);

    let (tip, base) = tip_and_base(&svg("c2:15-40")); // the minus gene starts at 11, before the left edge
    assert!((tip.0 - LEFT).abs() < 0.2 && base > tip.0);
}

#[test]
fn the_arrow_line_runs_along_the_visible_gene_below_the_gene_blocks() {
    let s = svg("c1:1-200");

    let x1: f64 = attr(&s, "class=\"gene-arrow\"", "x1").parse().unwrap();
    let x2: f64 = attr(&s, "class=\"gene-arrow\"", "x2").parse().unwrap();
    let y: f64 = attr(&s, "class=\"gene-arrow\"", "y1").parse().unwrap();
    assert!((x1 - (LEFT + 20.0 / 200.0 * PLOT)).abs() < 0.2 && (x2 - (LEFT + 110.0 / 200.0 * PLOT)).abs() < 0.2);
    let cds_top: f64 = attr(&s, "class=\"cds\"", "y").parse().unwrap();
    let cds_height: f64 = attr(&s, "class=\"cds\"", "height").parse().unwrap();
    assert!(y > cds_top + cds_height, "the arrow line {y} is below the blocks that end at {}", cds_top + cds_height);
}

#[test]
fn chevrons_along_the_gene_point_in_the_direction_of_transcription() {
    let plus = svg("c1:1-200");
    let minus = svg("c2:1-40");

    assert!(count(&plus, "class=\"chevron\"") >= 4, "a long gene has several chevrons");
    assert_eq!(attr(&plus, "class=\"chevron\"", "data-strand"), "+");
    assert_eq!(attr(&minus, "class=\"chevron\"", "data-strand"), "-");
    // A chevron is M x0,y-3 L x1,y L x0,y+3: for plus the tip (x1) is right of the open end (x0), for minus it is left.
    let d = |s: &str| -> (f64, f64) {
        let d = attr(s, "class=\"chevron\"", " d"); // with the space: `data-strand="` also ends in `d="`
        let nums: Vec<f64> = d.replace(['M', 'L', ','], " ").split_whitespace().map(|n| n.parse().unwrap()).collect();
        (nums[0], nums[2])
    };
    let (open, tip) = d(&plus);
    assert!(tip > open);
    let (open, tip) = d(&minus);
    assert!(tip < open);
}

#[test]
fn one_arrow_and_one_set_of_chevrons_per_transcript() {
    let s = svg("c1:1-200");

    assert_eq!(count(&s, "class=\"arrowhead\""), 1);
    assert_eq!(count(&s, "class=\"gene-arrow\""), 1);
    assert_eq!(count(&s, "class=\"intron\""), 1); // the existing intron line is unchanged
}

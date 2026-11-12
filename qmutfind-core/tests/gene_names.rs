//! Gene names and products from the GFF3 go into the annotated VCF (GENENAME, PRODUCT) and into the locus page.

use qmutfind_core::{annotate_text, read_report_variants, report_page};

/// A 200-base reference of repeated ACGT.
fn fasta() -> Vec<u8> {
    format!(">c1\n{}\n", "ACGT".repeat(50)).into_bytes()
}

const VCF: &str = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nc1\t30\t.\tG\tA\t.\tPASS\tSAMPLE=s\n";

fn gff(gene_attrs: &str, mrna_attrs: &str, cds_attrs: &str) -> String {
    format!(
        "##gff-version 3\n\
c1\tsrc\tgene\t21\t110\t.\t+\t.\tID=gene1;locus_tag=G1{gene_attrs}\n\
c1\tsrc\tmRNA\t21\t110\t.\t+\t.\tID=tx1;Parent=gene1{mrna_attrs}\n\
c1\tsrc\texon\t21\t110\t.\t+\t.\tParent=tx1\n\
c1\tsrc\tCDS\t21\t110\t.\t+\t0\tParent=tx1{cds_attrs}\n"
    )
}

fn info(gff_text: &str) -> String {
    let out = annotate_text(VCF, &fasta(), Some(gff_text.as_bytes())).unwrap();
    out.lines().last().unwrap().split('\t').nth(7).unwrap().to_string()
}

fn tag<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';').find_map(|kv| kv.strip_prefix(key)?.strip_prefix('='))
}

#[test]
fn a_gene_attribute_on_the_gene_line_gives_the_gene_name() {
    let i = info(&gff(";Name=hrmA;gene=hrmA", "", ""));

    assert_eq!(tag(&i, "GENE"), Some("G1"));
    assert_eq!(tag(&i, "GENENAME"), Some("hrmA"));
}

#[test]
fn a_name_that_differs_from_the_locus_tag_is_the_gene_name_when_there_is_no_gene_attribute() {
    assert_eq!(tag(&info(&gff(";Name=ssk1", "", "")), "GENENAME"), Some("ssk1"));
}

#[test]
fn a_name_equal_to_the_locus_tag_is_not_a_gene_name() {
    assert_eq!(tag(&info(&gff(";Name=G1", "", "")), "GENENAME"), None);
    assert_eq!(tag(&info(&gff("", "", "")), "GENENAME"), None);
}

#[test]
fn the_gene_attribute_wins_over_the_name_and_the_mrna_can_supply_the_name_when_the_gene_line_has_none() {
    assert_eq!(tag(&info(&gff(";Name=other;gene=hrmA", "", "")), "GENENAME"), Some("hrmA"));
    assert_eq!(tag(&info(&gff("", ";gene=fromMrna", "")), "GENENAME"), Some("fromMrna"));
}

#[test]
fn the_product_comes_from_the_mrna_or_else_from_the_cds() {
    assert_eq!(tag(&info(&gff("", ";product=response%20regulator", "")), "PRODUCT"), Some("response%20regulator"));
    assert_eq!(tag(&info(&gff("", "", ";product=kinase")), "PRODUCT"), Some("kinase"));
    assert_eq!(tag(&info(&gff("", ";product=from mrna", ";product=from cds")), "PRODUCT"), Some("from%20mrna"));
    assert_eq!(tag(&info(&gff("", "", "")), "PRODUCT"), None);
}

#[test]
fn characters_that_would_break_the_info_field_are_percent_encoded() {
    // GFF3 itself writes these as %3B, %3D and %2C; the VCF value must keep them encoded, and a space is encoded too.
    let i = info(&gff("", ";product=a%3Bb%3Dc%2Cd e%25f", ""));

    let product = tag(&i, "PRODUCT").unwrap();
    assert!(!product.contains(';') && !product.contains(' ') && !product.contains(','), "{product}");
    assert_eq!(i.matches('\t').count(), 0);
    assert_eq!(product, "a%3Bb%3Dc%2Cd%20e%25f");
}

#[test]
fn the_locus_page_table_and_drawing_use_the_gene_name() {
    let g = gff(";Name=hrmA;gene=hrmA", ";product=protein hrmA", "");
    let vcf = annotate_text(VCF, &fasta(), Some(g.as_bytes())).unwrap();

    let variants = read_report_variants(vcf.as_bytes()).unwrap();
    assert_eq!(variants[0].gene_name.as_deref(), Some("hrmA"));
    assert_eq!(variants[0].product.as_deref(), Some("protein hrmA"));

    let page = report_page(&vcf, Some(g.as_bytes()), "c1:1-200", "t").unwrap();
    assert!(page.contains("hrmA") && page.contains("protein hrmA"), "gene name and product are on the page");
    assert!(page.contains(">hrmA (G1)<"), "the gene model is labelled with the name and the ID");
}

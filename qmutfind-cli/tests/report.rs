//! End-to-end test of `qmutfind report`.

use std::fs;
use std::path::Path;
use std::process::Command;

const GFF: &str = "##gff-version 3
c1\tsrc\tgene\t21\t110\t.\t+\t.\tID=gene1;locus_tag=G1
c1\tsrc\tmRNA\t21\t110\t.\t+\t.\tID=tx1;Parent=gene1
c1\tsrc\texon\t21\t50\t.\t+\t.\tParent=tx1
c1\tsrc\texon\t81\t110\t.\t+\t.\tParent=tx1
c1\tsrc\tCDS\t24\t50\t.\t+\t0\tParent=tx1
c1\tsrc\tCDS\t81\t104\t.\t+\t0\tParent=tx1
";

const VCF: &str = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
c1\t28\t.\tC\tA\t.\tPASS\tEFFECT=missense;GENE=G1;TRANSCRIPT=tx1;AA=A2D
c1\t60\t.\tA\tC\t.\tPASS\tEFFECT=intron_derived;GENE=G1;TRANSCRIPT=tx1
c1\t500\t.\tT\tA\t.\tPASS\tEFFECT=intergenic
";

fn setup() -> tempfile::TempDir {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("in.vcf"), VCF).unwrap();
    fs::write(dir.path().join("genes.gff3"), GFF).unwrap();
    dir
}

fn run(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qmutfind")).arg("report").args(args).current_dir(dir).output().unwrap()
}

#[test]
fn report_writes_a_page_with_the_gene_the_variants_in_the_region_and_the_title() {
    let dir = setup();

    let out = run(dir.path(), &["--vcf", "in.vcf", "--gff", "genes.gff3", "--region", "c1:1-200", "--title", "My locus", "--out", "r.html"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let html = fs::read_to_string(dir.path().join("r.html")).unwrap();
    assert!(html.contains("<title>My locus</title>"));
    assert_eq!(html.matches("class=\"variant\"").count(), 2); // 28 and 60; 500 is outside the region
    assert_eq!(html.matches("class=\"cds\"").count(), 2);
}

#[test]
fn a_bad_region_is_reported_with_the_expected_form() {
    let dir = setup();

    let out = run(dir.path(), &["--vcf", "in.vcf", "--gff", "genes.gff3", "--region", "c1-200", "--out", "r.html"]);

    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("contig:start-end"));
}

#[test]
fn without_a_gff_the_page_still_shows_the_variants_and_no_gene_models() {
    let dir = setup();

    let out = run(dir.path(), &["--vcf", "in.vcf", "--region", "c1:1-200", "--out", "r.html"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let html = fs::read_to_string(dir.path().join("r.html")).unwrap();
    assert_eq!(html.matches("class=\"variant\"").count(), 2);
    assert_eq!(html.matches("class=\"cds\"").count(), 0);
}

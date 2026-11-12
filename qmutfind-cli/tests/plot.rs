//! End-to-end test of `qmutfind plot`.

use std::fs;
use std::path::Path;
use std::process::Command;

const VCF: &str = "##fileformat=VCFv4.2
##contig=<ID=c1,length=1000>
##contig=<ID=c2,length=2000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
c1\t100\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1.0000
c2\t500\t.\tCT\tC\t.\tPASS\tCLASS=SMALL_INDEL;AF=0.6000
";
const BARE: &str = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
c1\t100\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1.0000
";

fn run(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qmutfind")).arg("plot").args(args).current_dir(dir).output().unwrap()
}

#[test]
fn plot_writes_an_svg_with_one_marker_per_call() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("in.vcf"), VCF).unwrap();

    let out = run(dir.path(), &["--vcf", "in.vcf", "--out", "plot.svg"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let svg = fs::read_to_string(dir.path().join("plot.svg")).unwrap();
    assert!(svg.starts_with("<svg"));
    assert_eq!(svg.matches("class=\"call\"").count(), 2);
}

#[test]
fn plot_takes_contig_lengths_from_a_gzip_reference_when_the_vcf_has_none() {
    use std::io::Write;
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("in.vcf"), BARE).unwrap();
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(b">c1 first\nACGTACGTAC\n>c2\nACGT\n").unwrap();
    fs::write(dir.path().join("ref.fa.gz"), gz.finish().unwrap()).unwrap();

    let out = run(dir.path(), &["--vcf", "in.vcf", "--reference", "ref.fa.gz", "--out", "plot.svg"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(fs::read_to_string(dir.path().join("plot.svg")).unwrap().contains("data-contig-label=\"c2\""));
}

#[test]
fn plot_without_contig_lengths_fails_and_names_the_fix() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("in.vcf"), BARE).unwrap();

    let out = run(dir.path(), &["--vcf", "in.vcf", "--out", "plot.svg"]);

    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("reference"));
    assert!(!dir.path().join("plot.svg").exists());
}

#[test]
fn plot_orders_contigs_by_the_order_option_and_rejects_an_unknown_order() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("in.vcf"), VCF).unwrap(); // c1 is 1,000 bp and c2 is 2,000 bp

    let by_length = run(dir.path(), &["--vcf", "in.vcf", "--order", "length", "--out", "l.svg"]);
    assert!(by_length.status.success());
    let svg = fs::read_to_string(dir.path().join("l.svg")).unwrap();
    assert!(svg.find("data-contig-label=\"c2\"").unwrap() < svg.find("data-contig-label=\"c1\"").unwrap());

    let bad = run(dir.path(), &["--vcf", "in.vcf", "--order", "sideways", "--out", "b.svg"]);
    assert!(!bad.status.success());
}

fn depth_plot(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qmutfind")).arg("depth-plot").args(args).current_dir(dir).output().unwrap()
}

#[test]
fn depth_plot_reads_a_gzip_mosdepth_style_bed_and_draws_the_variants() {
    use std::io::Write;
    let dir = tempfile::tempdir().unwrap();
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(b"c1\t0\t500\t12.5\nc1\t500\t1000\t13.0\nc2\t0\t2000\t11.0\n").unwrap();
    fs::write(dir.path().join("sample.regions.bed.gz"), gz.finish().unwrap()).unwrap();
    fs::write(dir.path().join("in.vcf"), VCF).unwrap();

    let out = depth_plot(dir.path(), &["--bed", "sample.regions.bed.gz", "--vcf", "in.vcf", "--order", "name", "--out", "d.svg"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let svg = fs::read_to_string(dir.path().join("d.svg")).unwrap();
    assert!(svg.contains("class=\"depth\"") && svg.matches("class=\"vtick\"").count() == 2);
}

#[test]
fn depth_plot_reports_a_missing_bed_and_a_bad_row() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("bad.bed"), "c1\t0\tx\t5\n").unwrap();

    assert!(!depth_plot(dir.path(), &["--bed", "missing.bed", "--out", "d.svg"]).status.success());
    let out = depth_plot(dir.path(), &["--bed", "bad.bed", "--out", "d.svg"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("line 1"));
}

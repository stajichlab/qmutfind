//! End-to-end test of `qmutfind annotate`: VCF, GFF and reference in, annotated VCF out.
//! The gene model is the hand-checked fixture of qmutfind-core/tests/effects.rs (contig c1 only).

use std::fs;
use std::io::Write;
use std::path::Path;
use std::process::Command;

fn c1() -> String {
    format!(
        "{}CAC{}GT{}AG{}CCCCCC{}",
        "T".repeat(20),
        "ATGGCTGAACGTAAATTAGATCCCTGG",
        "A".repeat(26),
        "GGCAAACTGGACCATGAATGCTAA",
        "T".repeat(90)
    )
}

const GFF: &str = "##gff-version 3
c1\tsrc\tgene\t21\t110\t.\t+\t.\tID=gene1;locus_tag=G1
c1\tsrc\tmRNA\t21\t110\t.\t+\t.\tID=tx1;Parent=gene1
c1\tsrc\texon\t21\t50\t.\t+\t.\tParent=tx1
c1\tsrc\texon\t81\t110\t.\t+\t.\tParent=tx1
c1\tsrc\tCDS\t24\t50\t.\t+\t0\tParent=tx1
c1\tsrc\tCDS\t81\t104\t.\t+\t0\tParent=tx1
";

const VCF: &str = "##fileformat=VCFv4.2
##contig=<ID=c1,length=200>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
c1\t28\t.\tC\tA\t.\tPASS\tDP=20
c1\t60\t.\tA\tC\t.\tPASS\tDP=20
c1\t150\t.\tT\tA\t.\tPASS\tDP=20
";

fn setup() -> tempfile::TempDir {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("ref.fa"), format!(">c1\n{}\n", c1())).unwrap();
    fs::write(dir.path().join("genes.gff3"), GFF).unwrap();
    fs::write(dir.path().join("in.vcf"), VCF).unwrap();
    dir
}

fn run(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qmutfind")).arg("annotate").args(args).current_dir(dir).output().unwrap()
}

fn records(dir: &Path) -> Vec<String> {
    fs::read_to_string(dir.join("out.vcf")).unwrap().lines().filter(|l| !l.starts_with('#')).map(str::to_string).collect()
}

#[test]
fn annotate_adds_the_effect_tags_from_the_gff() {
    let dir = setup();

    let out = run(dir.path(), &["--vcf", "in.vcf", "--gff", "genes.gff3", "--reference", "ref.fa", "--out", "out.vcf"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert_eq!(
        records(dir.path()),
        vec![
            "c1\t28\t.\tC\tA\t.\tPASS\tDP=20;EFFECT=missense;GENE=G1;TRANSCRIPT=tx1;AA=A2D;NTX=1",
            "c1\t60\t.\tA\tC\t.\tPASS\tDP=20;EFFECT=intron_derived;GENE=G1;TRANSCRIPT=tx1;NTX=1",
            "c1\t150\t.\tT\tA\t.\tPASS\tDP=20;EFFECT=intergenic",
        ]
    );
}

#[test]
fn without_a_gff_every_record_is_unannotated_and_a_warning_is_printed() {
    let dir = setup();

    let out = run(dir.path(), &["--vcf", "in.vcf", "--reference", "ref.fa", "--out", "out.vcf"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(String::from_utf8_lossy(&out.stderr).contains("no GFF"));
    assert!(records(dir.path()).iter().all(|r| r.ends_with("DP=20;EFFECT=unannotated")));
}

#[test]
fn annotate_reads_gzip_vcf_and_gff() {
    let dir = setup();
    for (plain, gz) in [("in.vcf", "in.vcf.gz"), ("genes.gff3", "genes.gff3.gz")] {
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(&fs::read(dir.path().join(plain)).unwrap()).unwrap();
        fs::write(dir.path().join(gz), enc.finish().unwrap()).unwrap();
    }

    let out = run(dir.path(), &["--vcf", "in.vcf.gz", "--gff", "genes.gff3.gz", "--reference", "ref.fa", "--out", "out.vcf"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(records(dir.path())[0].contains("EFFECT=missense"));
}

#[test]
fn a_missing_required_option_is_named_in_the_error() {
    let dir = setup();

    let out = run(dir.path(), &["--gff", "genes.gff3", "--reference", "ref.fa", "--out", "out.vcf"]);

    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("missing required option --vcf"));
}

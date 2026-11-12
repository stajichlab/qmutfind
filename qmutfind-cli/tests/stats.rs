//! End-to-end test of `qmutfind stats`.

use std::fs;
use std::path::Path;
use std::process::Command;

const VCF: &str = "##fileformat=VCFv4.2
##contig=<ID=c1,length=100000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
c1\t100\t.\tA\tG\t.\tPASS\tCLASS=SNP;AF=1.0000
c1\t180\t.\tC\tT\t.\tPASS\tCLASS=SNP;AF=0.9000
c1\t5000\t.\tG\tT\t.\tPASS\tCLASS=SNP;AF=0.9500
";

fn run(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qmutfind")).arg("stats").args(args).current_dir(dir).output().unwrap()
}

#[test]
fn stats_writes_a_json_file_and_a_tsv_file_with_the_prefix() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("in.vcf"), VCF).unwrap();

    let out = run(dir.path(), &["--vcf", "in.vcf", "--out-prefix", "run1"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let json: serde_json::Value = serde_json::from_str(&fs::read_to_string(dir.path().join("run1.stats.json")).unwrap()).unwrap();
    assert_eq!(json["total"], 3);
    assert_eq!(json["transitions"], 2);
    assert!(fs::read_to_string(dir.path().join("run1.stats.tsv")).unwrap().starts_with("metric\tkey\tvalue\n"));
}

#[test]
fn the_clustering_window_defaults_to_100_bases_and_window_overrides_it() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("in.vcf"), VCF).unwrap(); // 100 and 180 are 80 apart

    assert!(run(dir.path(), &["--vcf", "in.vcf", "--out-prefix", "a"]).status.success());
    assert!(run(dir.path(), &["--vcf", "in.vcf", "--out-prefix", "b", "--window", "50"]).status.success());

    let clustered = |prefix: &str| -> serde_json::Value {
        serde_json::from_str::<serde_json::Value>(&fs::read_to_string(dir.path().join(format!("{prefix}.stats.json"))).unwrap()).unwrap()["clustered"].clone()
    };
    assert_eq!(clustered("a"), 2); // 80 bp apart is within 100
    assert_eq!(clustered("b"), 0); // but not within 50
}

#[test]
fn a_missing_vcf_option_is_named_in_the_error() {
    let dir = tempfile::tempdir().unwrap();

    let out = run(dir.path(), &["--out-prefix", "x"]);

    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("missing required option --vcf"));
}

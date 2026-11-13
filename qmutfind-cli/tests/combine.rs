//! End-to-end test of `qmutfind combine`.

use std::fs;
use std::path::Path;
use std::process::Command;

fn vcf(records: &[(u64, &str, &str, u64, &str)]) -> String {
    let mut text = String::from("##fileformat=VCFv4.2\n##contig=<ID=chrA,length=3000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (pos, reference, alt, dp, af) in records {
        text.push_str(&format!("chrA\t{pos}\t.\t{reference}\t{alt}\t50.0\tPASS\tSAMPLE=s;CLASS=SNP;DP={dp};AO=1;AF={af}\n"));
    }
    text
}

fn bed(depth: f64) -> String {
    (0..3).map(|i| format!("chrA\t{}\t{}\t{depth}\n", i * 1000, (i + 1) * 1000)).collect()
}

fn setup(dir: &Path, depth: f64) {
    fs::write(dir.join("parent.vcf"), vcf(&[(500, "A", "G", 30, "1.0000")])).unwrap();
    fs::write(dir.join("evolved.vcf"), vcf(&[(500, "A", "G", 28, "1.0000"), (1500, "C", "T", 25, "0.9600")])).unwrap();
    fs::write(dir.join("parent.bed"), bed(depth)).unwrap();
    fs::write(dir.join("evolved.bed"), bed(depth)).unwrap();
}

fn run(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qmutfind")).arg("combine").args(args).current_dir(dir).output().unwrap()
}

#[test]
fn combine_writes_the_table_and_the_summary_with_the_parent_statuses() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), 30.0);

    let out = run(dir.path(), &["--strain", "parent=parent.vcf,parent.bed", "--strain", "evolved=evolved.vcf,evolved.bed", "--parent", "parent", "--out-prefix", "run1"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let table = fs::read_to_string(dir.path().join("run1.combined.tsv")).unwrap();
    let lines: Vec<&str> = table.lines().collect();
    assert_eq!(lines[0], "CHROM\tPOS\tREF\tALT\tCLASS\tstrains\tparent\tevolved\tparent.status\tevolved.status");
    assert_eq!(lines[1], "chrA\t500\tA\tG\tSNP\t2\t1.00:30\t1.00:28\tparent\tinherited");
    assert_eq!(lines[2], "chrA\t1500\tC\tT\tSNP\t1\tref\t0.96:25\tparent\tnew");
    let summary = fs::read_to_string(dir.path().join("run1.summary.tsv")).unwrap();
    assert!(summary.contains("evolved\t2\t1\t1\t0\t0"), "{summary}");
}

#[test]
fn a_strain_without_a_bed_has_unknown_cells() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), 30.0);

    let out = run(dir.path(), &["--strain", "parent=parent.vcf", "--strain", "evolved=evolved.vcf", "--parent", "parent", "--out-prefix", "r"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let table = fs::read_to_string(dir.path().join("r.combined.tsv")).unwrap();
    assert!(table.lines().nth(2).unwrap().starts_with("chrA\t1500\tC\tT\tSNP\t1\tunknown\t0.96:25"), "{table}");
    assert!(table.lines().nth(2).unwrap().ends_with("parent\tunknown"), "{table}");
}

#[test]
fn min_depth_defaults_to_10_and_the_option_changes_it() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), 5.0); // window depth 5
    let args = |extra: &[&'static str]| {
        let mut a = vec!["--strain", "parent=parent.vcf,parent.bed", "--strain", "evolved=evolved.vcf,evolved.bed", "--out-prefix", "m"];
        a.extend_from_slice(extra);
        a
    };

    assert!(run(dir.path(), &args(&[])).status.success());
    assert!(fs::read_to_string(dir.path().join("m.combined.tsv")).unwrap().lines().nth(2).unwrap().contains("\tlow\t"));
    assert!(run(dir.path(), &args(&["--min-depth", "5"])).status.success());
    assert!(fs::read_to_string(dir.path().join("m.combined.tsv")).unwrap().lines().nth(2).unwrap().contains("\tref\t"));
}

#[test]
fn bad_input_fails_with_a_message_and_writes_nothing() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), 30.0);
    let cases: [(&[&str], &str); 4] = [
        (&["--strain", "nopath", "--out-prefix", "x"], "NAME=VCF"),
        (&["--strain", "a=parent.vcf", "--strain", "b=evolved.vcf", "--parent", "zzz", "--out-prefix", "x"], "parent"),
        (&["--strain", "a=missing.vcf", "--out-prefix", "x"], "missing.vcf"),
        (&["--out-prefix", "x"], "--strain"),
    ];
    for (args, fragment) in cases {
        let out = run(dir.path(), args);
        assert!(!out.status.success(), "{args:?}");
        assert!(String::from_utf8_lossy(&out.stderr).contains(fragment), "{args:?}: {}", String::from_utf8_lossy(&out.stderr));
    }
    assert!(!dir.path().join("x.combined.tsv").exists());
}

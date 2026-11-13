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
    // The merged VCF: one record per distinct variant, with the samples that call it and where it is new.
    let merged = fs::read_to_string(dir.path().join("run1.combined.vcf")).unwrap();
    let records: Vec<&str> = merged.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(records.len(), 2, "{merged}");
    assert!(records[0].contains("NS=2") && records[0].contains("SAMPLES=parent,evolved"), "{}", records[0]);
    assert!(records[1].contains("NS=1") && records[1].contains("NEW=evolved"), "{}", records[1]);
    assert!(merged.contains("##contig=<ID=chrA,length=3000>") && merged.contains("##qmutfind_parent=parent"));
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
fn min_depth_defaults_to_5_and_the_option_changes_it() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), 4.0); // window depth 4
    let args = |extra: &[&'static str]| {
        let mut a = vec!["--strain", "parent=parent.vcf,parent.bed", "--strain", "evolved=evolved.vcf,evolved.bed", "--out-prefix", "m"];
        a.extend_from_slice(extra);
        a
    };

    assert!(run(dir.path(), &args(&[])).status.success());
    assert!(fs::read_to_string(dir.path().join("m.combined.tsv")).unwrap().lines().nth(2).unwrap().contains("\tlow\t"));
    assert!(run(dir.path(), &args(&["--min-depth", "4"])).status.success());
    assert!(fs::read_to_string(dir.path().join("m.combined.tsv")).unwrap().lines().nth(2).unwrap().contains("\tref\t"));
}

#[test]
fn a_third_file_after_the_bed_is_the_evidence_vcf_and_an_empty_bed_field_is_allowed() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), 30.0);
    // The parent has 6 of 30 reads for the evolved strain's SNP at 1500, below the calling thresholds.
    fs::write(dir.path().join("parent.evidence.vcf"), vcf(&[(1500, "C", "T", 30, "0.2000")])).unwrap();

    let out = run(dir.path(), &["--strain", "parent=parent.vcf,parent.bed,parent.evidence.vcf", "--strain", "evolved=evolved.vcf,,", "--parent", "parent", "--out-prefix", "e"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let table = fs::read_to_string(dir.path().join("e.combined.tsv")).unwrap();
    let row = table.lines().nth(2).unwrap();
    assert!(row.starts_with("chrA\t1500\tC\tT\tSNP\t1\tweak:0.20:30\t0.96:25"), "{row}");
    assert!(row.ends_with("parent\tunknown"), "{row}");
}

#[test]
fn a_warning_about_a_coverage_contig_is_printed_to_stderr() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), 30.0);
    fs::write(dir.path().join("wrong.bed"), "chrZ\t0\t1000\t30\n").unwrap();

    let out = run(dir.path(), &["--strain", "a=parent.vcf,wrong.bed", "--out-prefix", "w"]);

    assert!(out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("chrZ"), "{}", String::from_utf8_lossy(&out.stderr));
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

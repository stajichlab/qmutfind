//! End-to-end test of `qmutfind call`: files in, VCF out.

use std::fs;
use std::path::Path;
use std::process::Command;

struct Lcg(u64);

impl Lcg {
    fn base(&mut self) -> u8 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        b"ACGT"[((self.0 >> 33) % 4) as usize]
    }
}

fn reference() -> Vec<u8> {
    let mut rng = Lcg(42);
    (0..60_000).map(|_| rng.base()).collect()
}

fn write_reference(dir: &Path, seq: &[u8]) {
    fs::write(dir.join("ref.fa"), format!(">chrA test\n{}\n", String::from_utf8_lossy(seq))).unwrap();
}

/// Ten 150 bp reads that all carry `alt` at 0-based position 30_000.
fn write_snp_reads(dir: &Path, seq: &[u8], alt: u8) {
    let mut fastq = String::new();
    for i in 0..10 {
        let start = 29_900 + i * 10;
        let mut read = seq[start..start + 150].to_vec();
        read[30_000 - start] = alt;
        fastq.push_str(&format!("@r{i}\n{}\n+\n{}\n", String::from_utf8_lossy(&read), "I".repeat(150)));
    }
    fs::write(dir.join("reads.fq"), fastq).unwrap();
}

fn records(vcf: &str) -> Vec<Vec<String>> {
    vcf.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect()
}

#[test]
fn call_writes_a_vcf_with_the_snp_carried_by_all_reads() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference();
    write_reference(dir.path(), &seq);
    let alt = if seq[30_000] == b'G' { b'T' } else { b'G' };
    write_snp_reads(dir.path(), &seq, alt);

    let status = Command::new(env!("CARGO_BIN_EXE_qmutfind"))
        .args(["call", "--reference", "ref.fa", "--reads", "reads.fq", "--sample", "s1"])
        .args(["--preset", "short", "--out", "out.vcf"])
        .current_dir(dir.path())
        .status()
        .unwrap();

    assert!(status.success());
    let vcf = fs::read_to_string(dir.path().join("out.vcf")).unwrap();
    let recs = records(&vcf);
    assert_eq!(recs.len(), 1, "{vcf}");
    assert_eq!(recs[0][..5], ["chrA", "30001", ".", &(seq[30_000] as char).to_string(), &(alt as char).to_string()]);
    assert!(recs[0][7].starts_with("SAMPLE=s1;CLASS=SNP;DP=10;AO=10;AF=1.0000"), "{}", recs[0][7]);
}

fn run_call(dir: &Path, reads: &str) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qmutfind"))
        .args(["call", "--reference", "ref.fa", "--reads", reads, "--sample", "s1"])
        .args(["--preset", "short", "--out", "out.vcf"])
        .current_dir(dir)
        .output()
        .unwrap()
}

fn setup_with_snp_reads() -> (tempfile::TempDir, Vec<u8>, u8) {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference();
    write_reference(dir.path(), &seq);
    let alt = if seq[30_000] == b'G' { b'T' } else { b'G' };
    write_snp_reads(dir.path(), &seq, alt);
    (dir, seq, alt)
}

#[test]
fn call_reads_gzip_compressed_fastq() {
    use std::io::Write;
    let (dir, _seq, _alt) = setup_with_snp_reads();
    let plain = fs::read(dir.path().join("reads.fq")).unwrap();
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    enc.write_all(&plain).unwrap();
    fs::write(dir.path().join("reads.fq.gz"), enc.finish().unwrap()).unwrap();

    let out = run_call(dir.path(), "reads.fq.gz");

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let vcf = fs::read_to_string(dir.path().join("out.vcf")).unwrap();
    assert_eq!(records(&vcf).len(), 1, "{vcf}");
}

#[test]
fn call_reads_zstd_compressed_fastq() {
    let (dir, _seq, _alt) = setup_with_snp_reads();
    let plain = fs::read(dir.path().join("reads.fq")).unwrap();
    fs::write(dir.path().join("reads.fq.zst"), zstd::encode_all(&plain[..], 3).unwrap()).unwrap();

    let out = run_call(dir.path(), "reads.fq.zst");

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let vcf = fs::read_to_string(dir.path().join("out.vcf")).unwrap();
    assert_eq!(records(&vcf).len(), 1, "{vcf}");
}

#[test]
fn a_corrupt_zstd_file_fails_with_a_message_and_a_non_zero_exit_status() {
    let (dir, _seq, _alt) = setup_with_snp_reads();
    fs::write(dir.path().join("bad.fq.zst"), b"this is not zstd data").unwrap();

    let out = run_call(dir.path(), "bad.fq.zst");

    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("bad.fq.zst"));
}

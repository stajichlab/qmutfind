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

fn run_with(dir: &Path, extra: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qmutfind"))
        .args(["call", "--reference", "ref.fa", "--reads", "reads.fq", "--sample", "s1"])
        .args(["--preset", "short", "--out", "out.vcf"])
        .args(extra)
        .current_dir(dir)
        .output()
        .unwrap()
}

#[test]
fn an_unknown_aligner_or_caller_is_rejected_with_the_allowed_values() {
    let (dir, _seq, _alt) = setup_with_snp_reads();

    let aligner = run_with(dir.path(), &["--aligner", "bowtie"]);
    let caller = run_with(dir.path(), &["--caller", "gatk"]);

    assert!(!aligner.status.success());
    assert!(String::from_utf8_lossy(&aligner.stderr).contains("--aligner must be rammap or bwa-mem2"));
    assert!(!caller.status.success());
    assert!(String::from_utf8_lossy(&caller.stderr).contains("--caller must be builtin or freebayes"));
}

#[test]
fn freebayes_with_the_rammap_aligner_is_rejected_because_it_needs_a_bam_from_an_external_aligner() {
    let (dir, _seq, _alt) = setup_with_snp_reads();

    let out = run_with(dir.path(), &["--aligner", "rammap", "--caller", "freebayes"]);

    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("freebayes needs an external aligner"));
    assert!(!dir.path().join("out.vcf").exists());
}

/// The external-tool tests need bwa-mem2, samtools and freebayes on PATH.
/// Run them with `pixi run test`, which puts the pinned tools on PATH.
#[test]
fn bwa_mem2_with_the_builtin_caller_finds_the_same_snp_as_the_rammap_path() {
    let (dir, seq, alt) = setup_with_snp_reads();

    let out = run_with(dir.path(), &["--aligner", "bwa-mem2", "--caller", "builtin"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let vcf = fs::read_to_string(dir.path().join("out.vcf")).unwrap();
    let recs = records(&vcf);
    assert_eq!(recs.len(), 1, "{vcf}");
    assert_eq!(recs[0][..5], ["chrA", "30001", ".", &(seq[30_000] as char).to_string(), &(alt as char).to_string()]);
    assert!(recs[0][7].starts_with("SAMPLE=s1;CLASS=SNP;DP=10;AO=10;AF=1.0000"), "{}", recs[0][7]);
}

#[test]
fn the_bwa_mem2_backend_removes_its_temporary_files_when_it_finishes() {
    let (dir, _seq, _alt) = setup_with_snp_reads();
    let tmp = tempfile::tempdir().unwrap();

    let out = run_with(dir.path(), &["--aligner", "bwa-mem2", "--tmpdir", tmp.path().to_str().unwrap()]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert_eq!(fs::read_dir(tmp.path()).unwrap().count(), 0, "temporary files were left behind");
}

#[test]
fn a_missing_external_tool_is_named_in_the_error_with_a_hint() {
    let (dir, _seq, _alt) = setup_with_snp_reads();

    let out = Command::new(env!("CARGO_BIN_EXE_qmutfind"))
        .args(["call", "--reference", "ref.fa", "--reads", "reads.fq", "--sample", "s1"])
        .args(["--preset", "short", "--out", "out.vcf", "--aligner", "bwa-mem2"])
        .env("PATH", "/nonexistent")
        .current_dir(dir.path())
        .output()
        .unwrap();

    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("bwa-mem2 was not found on PATH"), "{stderr}");
    assert!(stderr.contains("pixi run"), "{stderr}");
}

#[test]
fn bwa_mem2_with_freebayes_finds_the_snp_and_labels_it_with_the_sample_name() {
    let (dir, seq, alt) = setup_with_snp_reads();

    let out = run_with(dir.path(), &["--aligner", "bwa-mem2", "--caller", "freebayes"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let vcf = fs::read_to_string(dir.path().join("out.vcf")).unwrap();
    let recs = records(&vcf);
    assert_eq!(recs.len(), 1, "{vcf}");
    assert_eq!(recs[0][..5], ["chrA", "30001", ".", &(seq[30_000] as char).to_string(), &(alt as char).to_string()]);
    assert!(recs[0][7].starts_with("SAMPLE=s1;CLASS=SNP;DP=10;AO=10;AF=1.0000"), "{}", recs[0][7]);
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// Ten read pairs from 300 bp fragments. Mate 1 is the first 100 bp, mate 2 is the reverse complement of
/// the last 100 bp. Only mate 2 covers position 30_000, and only there is the base `alt`.
fn write_pairs_with_the_snp_in_mate2(dir: &Path, seq: &[u8], alt: u8) {
    let (mut r1, mut r2) = (String::new(), String::new());
    for i in 1..=10 {
        let start = 29_700 + i * 8;
        let frag = &seq[start..start + 300];
        let mate1 = frag[..100].to_vec();
        let mut mate2_fwd = frag[200..300].to_vec();
        mate2_fwd[30_000 - (start + 200)] = alt;
        let mate2 = reverse_complement(&mate2_fwd);
        r1.push_str(&format!("@p{i}/1\n{}\n+\n{}\n", String::from_utf8_lossy(&mate1), "I".repeat(100)));
        r2.push_str(&format!("@p{i}/2\n{}\n+\n{}\n", String::from_utf8_lossy(&mate2), "I".repeat(100)));
    }
    fs::write(dir.join("reads.fq"), r1).unwrap();
    fs::write(dir.join("reads2.fq"), r2).unwrap();
}

fn assert_the_mate2_snp_is_called(dir: &Path, seq: &[u8], alt: u8, out: &std::process::Output) {
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let vcf = fs::read_to_string(dir.join("out.vcf")).unwrap();
    let recs = records(&vcf);
    assert_eq!(recs.len(), 1, "{vcf}");
    assert_eq!(recs[0][..5], ["chrA", "30001", ".", &(seq[30_000] as char).to_string(), &(alt as char).to_string()]);
    assert!(recs[0][7].starts_with("SAMPLE=s1;CLASS=SNP;DP=10;AO=10;AF=1.0000"), "{}", recs[0][7]);
}

#[test]
fn reads2_is_mapped_as_the_second_mate_by_bwa_mem2() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference();
    write_reference(dir.path(), &seq);
    let alt = if seq[30_000] == b'G' { b'T' } else { b'G' };
    write_pairs_with_the_snp_in_mate2(dir.path(), &seq, alt);

    let out = run_with(dir.path(), &["--aligner", "bwa-mem2", "--reads2", "reads2.fq"]);

    assert_the_mate2_snp_is_called(dir.path(), &seq, alt, &out);
}

#[test]
fn reads2_is_mapped_as_more_single_end_reads_by_the_rammap_path() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference();
    write_reference(dir.path(), &seq);
    let alt = if seq[30_000] == b'G' { b'T' } else { b'G' };
    write_pairs_with_the_snp_in_mate2(dir.path(), &seq, alt);

    let out = run_with(dir.path(), &["--reads2", "reads2.fq"]);

    assert_the_mate2_snp_is_called(dir.path(), &seq, alt, &out);
}

#[test]
fn bwa_mem2_accepts_zstd_compressed_reads() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference();
    write_reference(dir.path(), &seq);
    let alt = if seq[30_000] == b'G' { b'T' } else { b'G' };
    write_pairs_with_the_snp_in_mate2(dir.path(), &seq, alt);
    for name in ["reads.fq", "reads2.fq"] {
        let plain = fs::read(dir.path().join(name)).unwrap();
        fs::write(dir.path().join(format!("{name}.zst")), zstd::encode_all(&plain[..], 3).unwrap()).unwrap();
    }

    let out = Command::new(env!("CARGO_BIN_EXE_qmutfind"))
        .args(["call", "--reference", "ref.fa", "--reads", "reads.fq.zst", "--reads2", "reads2.fq.zst"])
        .args(["--sample", "s1", "--preset", "short", "--out", "out.vcf", "--aligner", "bwa-mem2"])
        .current_dir(dir.path())
        .output()
        .unwrap();

    assert_the_mate2_snp_is_called(dir.path(), &seq, alt, &out);
}

/// Ten reads covering position 30_000, of which the first `with_alt` carry `alt` there.
fn write_reads_with_alt_fraction(dir: &Path, seq: &[u8], alt: u8, with_alt: usize) {
    let mut fastq = String::new();
    for i in 0..10 {
        let start = 29_900 + i * 10;
        let mut read = seq[start..start + 150].to_vec();
        if i < with_alt {
            read[30_000 - start] = alt;
        }
        fastq.push_str(&format!("@r{i}\n{}\n+\n{}\n", String::from_utf8_lossy(&read), "I".repeat(150)));
    }
    fs::write(dir.join("reads.fq"), fastq).unwrap();
}

fn setup_with_alt_fraction(with_alt: usize) -> tempfile::TempDir {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference();
    write_reference(dir.path(), &seq);
    let alt = if seq[30_000] == b'G' { b'T' } else { b'G' };
    write_reads_with_alt_fraction(dir.path(), &seq, alt, with_alt);
    dir
}

fn called_records(dir: &Path) -> usize {
    records(&fs::read_to_string(dir.join("out.vcf")).unwrap()).len()
}

#[test]
fn the_default_snp_threshold_for_short_reads_is_0_8_and_min_snp_fraction_overrides_it() {
    let dir = setup_with_alt_fraction(7); // fraction 0.7

    assert!(run_with(dir.path(), &[]).status.success());
    assert_eq!(called_records(dir.path()), 0);

    assert!(run_with(dir.path(), &["--min-snp-fraction", "0.5"]).status.success());
    assert_eq!(called_records(dir.path()), 1);
}

#[test]
fn min_depth_overrides_the_default_depth_threshold() {
    let dir = setup_with_alt_fraction(10); // depth 10, fraction 1.0

    assert!(run_with(dir.path(), &[]).status.success());
    assert_eq!(called_records(dir.path()), 1);

    assert!(run_with(dir.path(), &["--min-depth", "12"]).status.success());
    assert_eq!(called_records(dir.path()), 0);
}

#[test]
fn a_bad_threshold_value_is_rejected_with_the_flag_name() {
    let dir = setup_with_alt_fraction(10);

    let high = run_with(dir.path(), &["--min-snp-fraction", "1.5"]);
    let text = run_with(dir.path(), &["--min-indel-fraction", "abc"]);
    let depth = run_with(dir.path(), &["--min-depth", "-3"]);

    assert!(!high.status.success());
    assert!(String::from_utf8_lossy(&high.stderr).contains("--min-snp-fraction must be a number from 0 to 1"));
    assert!(!text.status.success());
    assert!(String::from_utf8_lossy(&text.stderr).contains("--min-indel-fraction must be a number from 0 to 1"));
    assert!(!depth.status.success());
    assert!(String::from_utf8_lossy(&depth.stderr).contains("--min-depth must be a whole number"));
}


#[test]
fn min_mapq_above_the_largest_mapping_quality_removes_every_alignment_and_a_bad_value_is_rejected() {
    let dir = setup_with_alt_fraction(10);

    assert!(run_with(dir.path(), &[]).status.success());
    assert_eq!(called_records(dir.path()), 1);

    assert!(run_with(dir.path(), &["--min-mapq", "61"]).status.success());
    assert_eq!(called_records(dir.path()), 0);

    let bad = run_with(dir.path(), &["--min-mapq", "300"]);
    assert!(!bad.status.success());
    assert!(String::from_utf8_lossy(&bad.stderr).contains("--min-mapq must be a whole number from 0 to 255"));
}

#[test]
fn min_mapq_also_filters_the_bwa_mem2_path() {
    let dir = setup_with_alt_fraction(10);

    assert!(run_with(dir.path(), &["--aligner", "bwa-mem2", "--min-mapq", "61"]).status.success());

    assert_eq!(called_records(dir.path()), 0);
}

#[test]
fn min_mapq_is_passed_to_freebayes() {
    let dir = setup_with_alt_fraction(10);

    assert!(run_with(dir.path(), &["--aligner", "bwa-mem2", "--caller", "freebayes"]).status.success());
    assert_eq!(called_records(dir.path()), 1);

    assert!(run_with(dir.path(), &["--aligner", "bwa-mem2", "--caller", "freebayes", "--min-mapq", "61"]).status.success());
    assert_eq!(called_records(dir.path()), 0);
}

/// Ten reads of 150 bases that all carry a 2 bp deletion at 0-based positions 30_001 and 30_002.
fn write_deletion_reads(dir: &Path, seq: &[u8]) {
    let mut fastq = String::new();
    for i in 0..10 {
        let start = 29_900 + i * 10;
        let mut read = seq[start..start + 152].to_vec();
        let at = 30_001 - start;
        read.drain(at..at + 2);
        fastq.push_str(&format!("@d{i}\n{}\n+\n{}\n", String::from_utf8_lossy(&read), "I".repeat(150)));
    }
    fs::write(dir.join("reads.fq"), fastq).unwrap();
}

fn setup_with_deletion_reads() -> tempfile::TempDir {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference();
    write_reference(dir.path(), &seq);
    write_deletion_reads(dir.path(), &seq);
    dir
}

#[test]
fn min_indel_alt_reads_removes_an_indel_with_too_few_reads_on_every_backend() {
    let dir = setup_with_deletion_reads();
    let backends: [&[&str]; 3] = [
        &[],
        &["--aligner", "bwa-mem2"],
        &["--aligner", "bwa-mem2", "--caller", "freebayes"],
    ];

    for backend in backends {
        let mut with_floor_10: Vec<&str> = backend.to_vec();
        with_floor_10.extend(["--min-indel-alt-reads", "10"]);
        assert!(run_with(dir.path(), &with_floor_10).status.success(), "{backend:?}");
        assert_eq!(called_records(dir.path()), 1, "ten reads pass a floor of ten: {backend:?}");

        let mut with_floor_11: Vec<&str> = backend.to_vec();
        with_floor_11.extend(["--min-indel-alt-reads", "11"]);
        assert!(run_with(dir.path(), &with_floor_11).status.success(), "{backend:?}");
        assert_eq!(called_records(dir.path()), 0, "ten reads fail a floor of eleven: {backend:?}");
    }
}

#[test]
fn a_bad_min_indel_alt_reads_value_is_rejected_with_the_flag_name() {
    let dir = setup_with_deletion_reads();

    let out = run_with(dir.path(), &["--min-indel-alt-reads", "-2"]);

    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--min-indel-alt-reads must be a whole number"));
}

fn call_rows(dir: &Path, extra: &[&str]) -> Vec<Vec<String>> {
    let status = Command::new(env!("CARGO_BIN_EXE_qmutfind"))
        .args(["call", "--reference", "ref.fa", "--reads", "reads.fq", "--sample", "s1", "--preset", "short", "--out", "out.vcf"])
        .args(extra)
        .current_dir(dir)
        .status()
        .unwrap();
    assert!(status.success());
    records(&fs::read_to_string(dir.join("out.vcf")).unwrap())
}

#[test]
fn call_writes_a_qual_and_a_depth_ratio_for_each_call() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference();
    write_reference(dir.path(), &seq);
    write_snp_reads(dir.path(), &seq, if seq[30_000] == b'G' { b'T' } else { b'G' });

    let rows = call_rows(dir.path(), &[]);

    assert_eq!(rows.len(), 1);
    let qual: f64 = rows[0][5].parse().expect("QUAL is a number");
    assert!(qual > 20.0, "{qual}");
    assert!(rows[0][7].contains(";DPR="), "{}", rows[0][7]);
}

#[test]
fn min_qual_and_max_depth_ratio_drop_calls() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference();
    write_reference(dir.path(), &seq);
    write_snp_reads(dir.path(), &seq, if seq[30_000] == b'G' { b'T' } else { b'G' });

    assert_eq!(call_rows(dir.path(), &["--min-qual", "1000"]).len(), 0); // above the cap of 999
    assert_eq!(call_rows(dir.path(), &["--min-qual", "20"]).len(), 1);
    assert_eq!(call_rows(dir.path(), &["--max-depth-ratio", "0.01"]).len(), 0);
    assert_eq!(call_rows(dir.path(), &["--max-depth-ratio", "50"]).len(), 1);
}

#[test]
fn depth_out_writes_windowed_mean_depth_as_a_four_column_bed() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference();
    write_reference(dir.path(), &seq);
    write_snp_reads(dir.path(), &seq, if seq[30_000] == b'G' { b'T' } else { b'G' });

    call_rows(dir.path(), &["--depth-out", "depth.bed", "--depth-window", "10000"]);

    let bed = fs::read_to_string(dir.path().join("depth.bed")).unwrap();
    let rows: Vec<Vec<&str>> = bed.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows.len(), 6); // 60,000 bp in windows of 10,000
    assert_eq!(rows[0][..3], ["chrA", "0", "10000"]);
    let window_with_reads: f64 = rows[2][3].parse().unwrap(); // reads span 29,900 to 30,150
    assert!(window_with_reads > 0.0, "{window_with_reads}");
    assert_eq!(rows[0][3], "0.00");
}

#[test]
fn depth_out_is_rejected_with_freebayes_because_that_path_keeps_no_counts() {
    let dir = tempfile::tempdir().unwrap();
    let seq = reference();
    write_reference(dir.path(), &seq);
    write_snp_reads(dir.path(), &seq, b'G');

    let out = Command::new(env!("CARGO_BIN_EXE_qmutfind"))
        .args(["call", "--reference", "ref.fa", "--reads", "reads.fq", "--sample", "s1", "--preset", "short", "--out", "out.vcf"])
        .args(["--aligner", "bwa-mem2", "--caller", "freebayes", "--depth-out", "depth.bed"])
        .current_dir(dir.path())
        .output()
        .unwrap();

    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--depth-out"));
}

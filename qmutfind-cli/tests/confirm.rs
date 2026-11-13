//! `qmutfind confirm`: is an introduced cassette integrated where it should be, and intact?

use std::fs;
use std::path::Path;
use std::process::Command;

struct Lcg(u64);
impl Lcg {
    fn base(&mut self) -> u8 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        b"ACGT"[((self.0 >> 33) % 4) as usize]
    }
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 20
    }
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', o => o }).collect()
}

/// Writes ref.fa (chr, 120 kb), hph.fa (a 2 kb cassette), genes.gff3 (a gene at 50,001-51,500) and reads.fq (30x from the genome with the gene
/// replaced by the cassette).
fn setup(dir: &Path, reverse: bool) {
    let mut rng = Lcg(21);
    let genome: Vec<u8> = (0..120_000).map(|_| rng.base()).collect();
    let cassette: Vec<u8> = (0..2_000).map(|_| rng.base()).collect();
    fs::write(dir.join("ref.fa"), format!(">chr test genome\n{}\n", String::from_utf8_lossy(&genome))).unwrap();
    fs::write(dir.join("hph.fa"), format!(">hph\n{}\n", String::from_utf8_lossy(&cassette))).unwrap();
    fs::write(
        dir.join("genes.gff3"),
        "##gff-version 3\nchr\ts\tgene\t50001\t51500\t.\t+\t.\tID=gene-AFUA_9G00010;Name=AFUA_9G00010;locus_tag=AFUA_9G00010\n\
         chr\ts\tmRNA\t50001\t51500\t.\t+\t.\tID=rna-XM_1;Parent=gene-AFUA_9G00010\n",
    )
    .unwrap();
    let insert = if reverse { revcomp(&cassette) } else { cassette };
    let modified = [&genome[..50_000], &insert[..], &genome[51_500..]].concat();
    let mut fastq = String::new();
    let n = modified.len() * 30 / 150;
    for i in 0..n {
        let start = (rng.next() as usize) % (modified.len() - 150);
        let mut read = modified[start..start + 150].to_vec();
        if rng.next() & 1 == 1 {
            read = revcomp(&read);
        }
        fastq.push_str(&format!("@r{i}\n{}\n+\n{}\n", String::from_utf8_lossy(&read), "I".repeat(150)));
    }
    fs::write(dir.join("reads.fq"), fastq).unwrap();
}

fn confirm(dir: &Path, extra: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qmutfind"))
        .args(["confirm", "--reference", "ref.fa", "--construct", "hph.fa", "--reads", "reads.fq", "--preset", "short", "--out-prefix", "out"])
        .args(extra)
        .current_dir(dir)
        .output()
        .unwrap()
}

#[test]
fn confirm_reports_the_site_the_orientation_and_the_checks_and_writes_the_files() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), false);

    let out = confirm(dir.path(), &["--target", "chr:50001-51500"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let text = fs::read_to_string(dir.path().join("out.confirm.txt")).unwrap();
    assert!(text.contains("CONFIRMED") && !text.contains("NOT CONFIRMED"), "{text}");
    assert!(text.contains("left flank 50000") && text.contains("right flank 51501") && text.contains("orientation forward"), "{text}");
    let json = fs::read_to_string(dir.path().join("out.confirm.json")).unwrap();
    assert!(json.contains("\"verdict\":\"confirmed\""), "{json}");
    let junctions = fs::read_to_string(dir.path().join("out.junctions.tsv")).unwrap();
    assert!(junctions.lines().next().unwrap().starts_with("contig_a\tpos_a"), "{junctions}");
    assert!(junctions.lines().filter(|l| l.contains("hph")).count() >= 2, "{junctions}");
    let depth = fs::read_to_string(dir.path().join("out.cassette_depth.tsv")).unwrap();
    assert_eq!(depth.lines().count(), 1 + 2_000, "one line per cassette base and a header");
    assert!(String::from_utf8_lossy(&out.stdout).contains("CONFIRMED"), "the report is also printed");
}

#[test]
fn confirm_finds_a_gene_by_name_in_the_gff_and_a_reverse_cassette() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), true);

    let out = confirm(dir.path(), &["--gff", "genes.gff3", "--gene", "AFUA_9G00010"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let text = fs::read_to_string(dir.path().join("out.confirm.txt")).unwrap();
    assert!(text.contains("target chr:50001-51500"), "{text}");
    assert!(text.contains("orientation reverse") && text.contains("CONFIRMED"), "{text}");
}

#[test]
fn confirm_without_a_target_still_reports_the_site() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), false);

    let out = confirm(dir.path(), &[]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let text = fs::read_to_string(dir.path().join("out.confirm.txt")).unwrap();
    assert!(text.contains("left flank 50000"), "{text}");
    assert!(text.contains("[n/a] at_target"), "{text}");
}

#[test]
fn confirm_exits_with_a_status_that_says_whether_the_allele_was_confirmed() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), false);
    // The target is somewhere else: not confirmed. The command still writes its report, and exits with status 3.
    let out = confirm(dir.path(), &["--target", "chr:80001-81500"]);
    assert_eq!(out.status.code(), Some(3), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(fs::read_to_string(dir.path().join("out.confirm.txt")).unwrap().contains("NOT CONFIRMED"));
    // With --no-fail-status the exit status is 0 whatever the verdict.
    let out = confirm(dir.path(), &["--target", "chr:80001-81500", "--always-exit-zero"]);
    assert_eq!(out.status.code(), Some(0));
}

#[test]
fn confirm_refuses_bad_input_with_a_reason() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), false);
    let name_clash = fs::read_to_string(dir.path().join("hph.fa")).unwrap().replace(">hph", ">chr");
    fs::write(dir.path().join("clash.fa"), name_clash).unwrap();
    let run = |args: &[&str]| {
        let out = Command::new(env!("CARGO_BIN_EXE_qmutfind")).arg("confirm").args(args).current_dir(dir.path()).output().unwrap();
        (out.status.success(), String::from_utf8_lossy(&out.stderr).to_string())
    };
    let base = ["--reference", "ref.fa", "--reads", "reads.fq", "--preset", "short", "--out-prefix", "o"];
    let (ok, err) = run(&[&base[..], &["--construct", "clash.fa"]].concat());
    assert!(!ok && err.contains("chr") && err.contains("same name"), "{err}");
    let (ok, err) = run(&[&base[..], &["--construct", "hph.fa", "--target", "chr:50"]].concat());
    assert!(!ok && err.contains("CONTIG:START-END"), "{err}");
    let (ok, err) = run(&[&base[..], &["--construct", "hph.fa", "--gene", "AFUA_9G00010"]].concat());
    assert!(!ok && err.contains("--gff"), "{err}");
    let (ok, err) = run(&[&base[..], &["--construct", "hph.fa", "--gff", "genes.gff3", "--gene", "NOPE"]].concat());
    assert!(!ok && err.contains("NOPE"), "{err}");
    let (ok, err) = run(&["--reference", "ref.fa", "--reads", "reads.fq", "--preset", "short", "--out-prefix", "o"]);
    assert!(!ok && err.contains("--construct"), "{err}");
    let (ok, err) = run(&[&base[..], &["--construct", "hph.fa", "--target", "nope:1-100"]].concat());
    assert!(!ok && err.contains("nope"), "{err}");
}

/// A donor: 300 bases of the genome on each side of the cassette.
fn write_donor(dir: &Path) {
    let genome = fs::read_to_string(dir.join("ref.fa")).unwrap().lines().nth(1).unwrap().to_string();
    let cassette = fs::read_to_string(dir.join("hph.fa")).unwrap().lines().nth(1).unwrap().to_string();
    fs::write(dir.join("donor.fa"), format!(">donor\n{}{}{}\n", &genome[49_700..50_000], cassette, &genome[51_500..51_800])).unwrap();
}

#[test]
fn confirm_cuts_homology_arms_and_takes_the_target_from_them() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), false);
    write_donor(dir.path());

    let out = Command::new(env!("CARGO_BIN_EXE_qmutfind"))
        .args(["confirm", "--reference", "ref.fa", "--construct", "donor.fa", "--reads", "reads.fq", "--preset", "short", "--out-prefix", "d"])
        .current_dir(dir.path()).output().unwrap();

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let text = fs::read_to_string(dir.path().join("d.confirm.txt")).unwrap();
    assert!(text.contains("CONFIRMED") && text.contains("inferred from the homology arms") && text.contains("arms were cut off"), "{text}");
    // This test genome shares 4 bases by chance between the start of the cassette and the genome after the arm, so the arm edge, the target start
    // and the left flank are each uncertain by up to 4 bases (see the sweep in qmutfind-core/tests/confirm.rs).
    let line = text.lines().find(|l| l.starts_with("  chr donor")).unwrap_or_else(|| panic!("{text}"));
    assert!((50_000..=50_004).any(|p| line.contains(&format!("left flank {p} "))) && line.contains("right flank 51501"), "{text}");
    assert!((50_001..=50_005).any(|p| text.contains(&format!("target chr:{p}-51500"))), "{text}");
}

#[test]
fn confirm_keep_arms_keeps_them_and_warns() {
    let dir = tempfile::tempdir().unwrap();
    setup(dir.path(), false);
    write_donor(dir.path());

    let out = Command::new(env!("CARGO_BIN_EXE_qmutfind"))
        .args(["confirm", "--reference", "ref.fa", "--construct", "donor.fa", "--reads", "reads.fq", "--preset", "short", "--out-prefix", "d", "--keep-arms", "--always-exit-zero"])
        .current_dir(dir.path()).output().unwrap();

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let text = fs::read_to_string(dir.path().join("d.confirm.txt")).unwrap();
    assert!(text.contains("Warning:") && text.contains("homology arms") && !text.contains("arms were cut off"), "{text}");
}

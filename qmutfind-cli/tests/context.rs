//! `qmutfind context`: the sequence around a variant, as text or HTML.

use std::fs;
use std::process::Command;

struct Lcg(u64);
impl Lcg {
    fn base(&mut self) -> u8 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        b"ACGT"[((self.0 >> 33) % 4) as usize]
    }
}

/// A 200-base contig with a plus-strand gene whose CDS is 41..=70 and 91..=126 (ATG at 41, TAA at 124).
fn files(dir: &std::path::Path) -> Vec<u8> {
    let mut rng = Lcg(3);
    let mut seq: Vec<u8> = (0..200).map(|_| rng.base()).collect();
    seq[40..43].copy_from_slice(b"ATG");
    seq[123..126].copy_from_slice(b"TAA");
    fs::write(dir.join("ref.fa"), format!(">c1\n{}\n", String::from_utf8_lossy(&seq))).unwrap();
    fs::write(
        dir.join("genes.gff3"),
        "##gff-version 3\nc1\ts\tgene\t41\t126\t.\t+\t.\tID=g1;locus_tag=G1\nc1\ts\tmRNA\t41\t126\t.\t+\t.\tID=t1;Parent=g1\n\
         c1\ts\texon\t41\t70\t.\t+\t.\tParent=t1\nc1\ts\texon\t91\t126\t.\t+\t.\tParent=t1\n\
         c1\ts\tCDS\t41\t70\t.\t+\t0\tParent=t1\nc1\ts\tCDS\t91\t126\t.\t+\t.\tParent=t1\n",
    )
    .unwrap();
    seq
}

fn run(dir: &std::path::Path, extra: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qmutfind")).args(["context", "--reference", "ref.fa", "--gff", "genes.gff3"]).args(extra).current_dir(dir).output().unwrap()
}

#[test]
fn context_prints_the_aligned_rows_of_a_variant() {
    let dir = tempfile::tempdir().unwrap();
    let seq = files(dir.path());
    let r = seq[49] as char;
    let alt = if r == 'G' { 'T' } else { 'G' };

    let out = run(dir.path(), &["--variant", &format!("c1:50:{r}:{alt}"), "--flank", "10"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let text = String::from_utf8(out.stdout).unwrap();
    let line = |label: &str| text.lines().find(|l| l.starts_with(label)).unwrap_or_else(|| panic!("no {label} row in\n{text}")).to_string();
    assert!(text.lines().next().unwrap().starts_with("c1:50 "), "{text}");
    assert_eq!(line("marker").trim_end().find('^'), Some(12 + 10), "{text}");
    assert_eq!(line("reference").len(), 12 + 21);
    assert_eq!(line("variant").len(), line("reference").len());
}

#[test]
fn context_writes_html_with_the_changed_column_marked_when_asked() {
    let dir = tempfile::tempdir().unwrap();
    let seq = files(dir.path());
    let r = seq[49] as char;
    let alt = if r == 'G' { 'T' } else { 'G' };

    let out = run(dir.path(), &["--variant", &format!("c1:50:{r}:{alt}"), "--html", "--out", "ctx.html"]);

    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let html = fs::read_to_string(dir.path().join("ctx.html")).unwrap();
    assert!(html.contains("<span class=\"ctx-var\">"), "{html}");
}

#[test]
fn context_without_a_gff_shows_the_forward_strand_and_says_so() {
    let dir = tempfile::tempdir().unwrap();
    let seq = files(dir.path());
    let r = seq[49] as char;
    let out = Command::new(env!("CARGO_BIN_EXE_qmutfind"))
        .args(["context", "--reference", "ref.fa", "--variant", &format!("c1:50:{r}:A"), "--flank", "5"])
        .current_dir(dir.path()).output().unwrap();
    assert!(out.status.success() || r == 'A', "{}", String::from_utf8_lossy(&out.stderr));
    if r != 'A' {
        assert!(String::from_utf8_lossy(&out.stdout).contains("not in a gene"));
    }
}

#[test]
fn a_bad_variant_text_or_a_wrong_reference_allele_is_an_error_with_a_reason() {
    let dir = tempfile::tempdir().unwrap();
    let seq = files(dir.path());
    let bad = run(dir.path(), &["--variant", "c1:50"]);
    assert!(!bad.status.success());
    assert!(String::from_utf8_lossy(&bad.stderr).contains("CONTIG:POS:REF:ALT"), "{}", String::from_utf8_lossy(&bad.stderr));
    let wrong = if seq[49] == b'A' { 'C' } else { 'A' };
    let out = run(dir.path(), &["--variant", &format!("c1:50:{wrong}:T")]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("reference allele"), "{}", String::from_utf8_lossy(&out.stderr));
}

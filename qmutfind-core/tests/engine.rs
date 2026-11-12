//! The whole pipeline over bytes: reference FASTA and FASTQ chunks in, VCF text out.
//! The browser build and the CLI both drive this type, so they give the same result.

mod common;

use common::{other_base, reference};
use qmutfind_core::{Engine, ReadKind};

fn fasta(refs: &[(String, Vec<u8>)]) -> Vec<u8> {
    let mut out = Vec::new();
    for (name, seq) in refs {
        out.extend_from_slice(format!(">{name} description\n").as_bytes());
        for line in seq.chunks(60) {
            out.extend_from_slice(line);
            out.push(b'\n');
        }
    }
    out
}

/// FASTQ text for `count` reads covering 0-based position 50_000, the first `with_alt` carrying `alt`.
fn fastq(refs: &[(String, Vec<u8>)], alt: u8, with_alt: usize, count: usize, prefix: &str) -> Vec<u8> {
    let mut out = String::new();
    for i in 0..count {
        let start = 49_900 + i * 10;
        let mut read = refs[0].1[start..start + 150].to_vec();
        if i < with_alt {
            read[50_000 - start] = alt;
        }
        out.push_str(&format!("@{prefix}{i}\n{}\n+\n{}\n", String::from_utf8_lossy(&read), "I".repeat(150)));
    }
    out.into_bytes()
}

fn records(vcf: &str) -> Vec<Vec<String>> {
    vcf.lines().filter(|l| !l.starts_with('#')).map(|l| l.split('\t').map(str::to_string).collect()).collect()
}

#[test]
fn reads_fed_in_small_chunks_give_the_snp_record_the_cli_gives() {
    let refs = reference();
    let alt = other_base(refs[0].1[50_000]);
    let mut engine = Engine::new(&fasta(&refs), ReadKind::Short, "s1").unwrap();

    for chunk in fastq(&refs, alt, 10, 10, "r").chunks(7) {
        engine.add_reads(chunk).unwrap();
    }
    engine.end_file().unwrap();
    let vcf = engine.vcf().unwrap();

    let recs = records(&vcf);
    assert_eq!(recs.len(), 1, "{vcf}");
    assert_eq!(recs[0][..5], ["chrA", "50001", ".", &(refs[0].1[50_000] as char).to_string(), &(alt as char).to_string()]);
    assert!(recs[0][7].starts_with("SAMPLE=s1;CLASS=SNP;DP=10;AO=10;AF=1.0000"), "{}", recs[0][7]);
    assert!(vcf.contains("##contig=<ID=chrA,length=120000>"), "{vcf}");
}

#[test]
fn reads_from_two_files_are_both_counted_and_the_read_count_is_reported() {
    let refs = reference();
    let alt = other_base(refs[0].1[50_000]);
    let mut engine = Engine::new(&fasta(&refs), ReadKind::Short, "s1").unwrap();

    engine.add_reads(&fastq(&refs, alt, 3, 3, "a")).unwrap();
    engine.end_file().unwrap();
    engine.add_reads(&fastq(&refs, alt, 4, 4, "b")).unwrap();
    engine.end_file().unwrap();

    assert_eq!(engine.read_count(), 7);
    let recs = records(&engine.vcf().unwrap());
    assert_eq!(recs.len(), 1);
    assert!(recs[0][7].contains("DP=7;AO=7"), "{}", recs[0][7]);
}

#[test]
fn set_thresholds_changes_what_is_called() {
    let refs = reference();
    let alt = other_base(refs[0].1[50_000]);
    let data = fastq(&refs, alt, 7, 10, "r"); // fraction 0.7

    let mut strict = Engine::new(&fasta(&refs), ReadKind::Short, "s1").unwrap(); // short default: SNP 0.8
    strict.add_reads(&data).unwrap();
    strict.end_file().unwrap();
    assert!(records(&strict.vcf().unwrap()).is_empty());

    let mut lenient = Engine::new(&fasta(&refs), ReadKind::Short, "s1").unwrap();
    lenient.set_thresholds(5, 0.5, 0.5);
    lenient.add_reads(&data).unwrap();
    lenient.end_file().unwrap();
    assert_eq!(records(&lenient.vcf().unwrap()).len(), 1);
}

#[test]
fn a_truncated_read_file_is_an_error_at_end_file() {
    let refs = reference();
    let mut engine = Engine::new(&fasta(&refs), ReadKind::Short, "s1").unwrap();

    engine.add_reads(b"@r1\nACGT\n+\nIIII\n@r2\nGGGG\n").unwrap();

    assert!(engine.end_file().is_err());
}

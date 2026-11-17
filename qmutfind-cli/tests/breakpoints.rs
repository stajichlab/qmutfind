//! `qmutfind breakpoints` and `qmutfind locus`: find an insertion in a simulated assembly, write the table and the sequence, and draw it.

use std::fs;
use std::process::Command;

struct Lcg(u64);
impl Lcg {
    fn base(&mut self) -> u8 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        b"ACGT"[((self.0 >> 33) % 4) as usize]
    }
}

fn dna(n: usize, seed: u64) -> Vec<u8> {
    let mut rng = Lcg(seed);
    (0..n).map(|_| rng.base()).collect()
}

fn text(seq: &[u8]) -> String {
    String::from_utf8_lossy(seq).into_owned()
}

/// A 60 kb reference with a FungiDB-style gene (a `protein_coding_gene` line) at 19,001-23,000 and two assemblies: one with 2,500 bases replaced by 1,500 other bases, one unchanged.
fn files(dir: &std::path::Path) -> Vec<u8> {
    let genome = dna(60_000, 1);
    let insert = dna(1_500, 2);
    fs::write(dir.join("ref.fa"), format!(">chr\n{}\n", text(&genome))).unwrap();
    fs::write(dir.join("mut.fasta"), format!(">ctg1\n{}{}{}\n", text(&genome[..20_000]), text(&insert), text(&genome[22_500..]))).unwrap();
    fs::write(dir.join("wt.fasta"), format!(">ctg1\n{}\n", text(&genome))).unwrap();
    fs::write(
        dir.join("genes.gff3"),
        "##gff-version 3\nchr\ts\tprotein_coding_gene\t19001\t23000\t.\t+\t.\tID=G1\nchr\ts\tmRNA\t19001\t23000\t.\t+\t.\tID=G1-T;Parent=G1\n\
         chr\ts\texon\t19001\t23000\t.\t+\t.\tParent=G1-T\nchr\ts\tCDS\t19101\t22900\t.\t+\t0\tParent=G1-T\n",
    )
    .unwrap();
    insert
}

fn cli(dir: &std::path::Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qmutfind")).args(args).current_dir(dir).output().unwrap()
}

#[test]
fn an_assembly_is_aligned_and_the_gap_the_table_the_sequence_and_the_drawing_are_written() {
    let dir = std::env::temp_dir().join(format!("qmutfind-bp-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let insert = files(&dir);

    let out = cli(&dir, &["breakpoints", "--reference", "ref.fa", "--gff", "genes.gff3", "--gene", "G1", "--assembly", "mutant=mut.fasta", "--assembly", "wt.fasta", "--out-prefix", "bp"]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(stdout.contains("mutant: replacement of 15") && stdout.contains("between chr:200"), "{stdout}");
    assert!(stdout.contains("wt: alignments present, no gap overlapping chr:19001-23000"), "{stdout}");

    let table = fs::read_to_string(dir.join("bp.breakpoint_summary.tsv")).unwrap();
    let rows: Vec<Vec<&str>> = table.lines().skip(1).map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows.len(), 1, "{table}");
    assert_eq!((rows[0][0], rows[0][3], rows[0][8], rows[0][9]), ("mutant", "chr", "replacement", "yes"));
    let fasta = fs::read_to_string(dir.join("bp.inserted_sequences.fasta")).unwrap();
    let seq = fasta.lines().nth(1).unwrap();
    assert!(text(&insert).contains(&seq[3..seq.len() - 3]), "the inserted sequence is the simulated one");
    // The gene is found as a protein_coding_gene, the insertion removes its coding sequence, and the drawing is written.
    let genes = fs::read_to_string(dir.join("bp.mutant.locus_genes.tsv")).unwrap();
    assert!(genes.contains("G1\t") && genes.contains("chr:19001-23000"), "{genes}");
    assert!(fs::read_to_string(dir.join("bp.mutant.locus.svg")).unwrap().starts_with("<svg"));

    // `locus` draws a row of the same table, with a link for the hit gene.
    let out = cli(&dir, &["locus", "--breakpoints", "bp.breakpoint_summary.tsv", "--gff", "genes.gff3", "--gene-link", "https://fungidb.org/gene/", "--out-prefix", "again"]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(fs::read_to_string(dir.join("again.locus.svg")).unwrap().contains("href=\"https://fungidb.org/gene/G1\""));
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_missing_locus_assembly_or_gene_is_an_error_that_says_what_to_give() {
    let dir = std::env::temp_dir().join(format!("qmutfind-bp-err-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    files(&dir);
    let reason = |args: &[&str]| String::from_utf8_lossy(&cli(&dir, args).stderr).into_owned();
    assert!(reason(&["breakpoints", "--reference", "ref.fa", "--out-prefix", "x"]).contains("--assembly"));
    assert!(reason(&["breakpoints", "--reference", "ref.fa", "--assembly", "mut.fasta", "--out-prefix", "x"]).contains("give the locus"));
    assert!(reason(&["breakpoints", "--reference", "ref.fa", "--assembly", "mut.fasta", "--gene", "G1", "--out-prefix", "x"]).contains("--gene needs --gff"));
    assert!(reason(&["breakpoints", "--reference", "ref.fa", "--assembly", "mut.fasta", "--locus", "nope:1-10", "--out-prefix", "x"]).contains("nope"));
    assert!(reason(&["breakpoints", "--reference", "ref.fa", "--assembly", "mut.fasta", "--locus", "chr:1-10", "--preset", "asm7", "--out-prefix", "x"]).contains("asm5, asm10 or asm20"));
    fs::remove_dir_all(&dir).unwrap();
}

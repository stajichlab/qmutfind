//! VCF output for calls. One record per sample per mutation, with the sample in INFO.

use qmutfind_core::{write_vcf, Call};

fn call(pos: u64, reference: &str, alt: &str, depth: u32, alt_count: u32) -> Call {
    Call { contig: "chrA".to_string(), pos, reference: reference.to_string(), alt: alt.to_string(), depth, alt_count }
}

fn render(calls: &[Call]) -> Vec<String> {
    let contigs = vec![("chrA".to_string(), 1000usize), ("chrB".to_string(), 500usize)];
    let mut out = Vec::new();
    write_vcf(&mut out, &contigs, "mutant1", calls, 50).unwrap();
    String::from_utf8(out).unwrap().lines().map(str::to_string).collect()
}

#[test]
fn the_header_names_the_format_the_contigs_and_every_info_field() {
    let lines = render(&[]);

    assert_eq!(lines[0], "##fileformat=VCFv4.2");
    assert!(lines.contains(&"##contig=<ID=chrA,length=1000>".to_string()));
    assert!(lines.contains(&"##contig=<ID=chrB,length=500>".to_string()));
    for id in ["SAMPLE", "CLASS", "DP", "AC", "AF"] {
        assert!(lines.iter().any(|l| l.starts_with(&format!("##INFO=<ID={id},"))), "missing INFO {id}");
    }
    assert_eq!(lines.last().unwrap(), "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
}

#[test]
fn each_call_is_one_record_with_sample_class_depth_and_allele_fraction() {
    let big_del = format!("C{}", "A".repeat(50));
    let lines = render(&[
        call(11, "A", "G", 10, 10),
        call(21, "CTT", "C", 8, 6),
        call(31, &big_del, "C", 4, 3),
    ]);

    let records: Vec<&String> = lines.iter().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(records.len(), 3);
    assert_eq!(records[0].as_str(), "chrA\t11\t.\tA\tG\t.\tPASS\tSAMPLE=mutant1;CLASS=SNP;DP=10;AC=10;AF=1.0000");
    assert_eq!(records[1].as_str(), "chrA\t21\t.\tCTT\tC\t.\tPASS\tSAMPLE=mutant1;CLASS=SMALL_INDEL;DP=8;AC=6;AF=0.7500");
    assert_eq!(records[2].as_str(), format!("chrA\t31\t.\t{big_del}\tC\t.\tPASS\tSAMPLE=mutant1;CLASS=LARGE_INDEL;DP=4;AC=3;AF=0.7500"));
}

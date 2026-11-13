//! The whole path for a simulated parent and child: reads -> caller -> VCF, coverage and evidence -> combine.
//! Each site has a known truth, so the status of every row is known.

mod common;

use common::{other_base, reference};
use qmutfind_core::{combine, parse_sam_line, write_vcf, Caller, CallerConfig, Cell, CombineInput, CombineOptions, Status};

/// One site: 0-based position, the reads that carry the alternate base, and the depth.
struct Site {
    pos: usize,
    alt_reads: usize,
    depth: usize,
}

struct Outputs {
    vcf: String,
    bed: String,
    evidence: String,
}

fn run_strain(sites: &[Site]) -> Outputs {
    let refs = reference();
    let config = CallerConfig { min_depth: 5, min_snp_fraction: 0.8, min_indel_fraction: 0.8, large_indel_min: 50 };
    let mut caller = Caller::new(refs.clone(), config);
    for site in sites {
        let alt = other_base(refs[0].1[site.pos]);
        for i in 0..site.depth {
            let start = site.pos - 75 + (i % 20);
            let mut seq = refs[0].1[start..start + 150].to_vec();
            if i < site.alt_reads {
                seq[site.pos - start] = alt;
            }
            let sam = format!("r{i}\t0\tchrA\t{}\t60\t150M\t*\t0\t0\t{}\t{}", start + 1, String::from_utf8(seq).unwrap(), "I".repeat(150));
            caller.add_sam(&parse_sam_line(&sam).unwrap().unwrap());
        }
    }
    let contigs: Vec<(String, usize)> = refs.iter().map(|(n, s)| (n.clone(), s.len())).collect();
    let mut vcf = Vec::new();
    write_vcf(&mut vcf, &contigs, "s", &caller.calls(), 50).unwrap();
    let mut evidence = Vec::new();
    write_vcf(&mut evidence, &contigs, "s", &caller.evidence_calls(), 50).unwrap();
    Outputs { vcf: String::from_utf8(vcf).unwrap(), bed: caller.depth_bed(100), evidence: String::from_utf8(evidence).unwrap() }
}

fn site(pos: usize, alt_reads: usize, depth: usize) -> Site {
    Site { pos, alt_reads, depth }
}

#[test]
fn a_simulated_parent_and_child_give_the_known_status_at_every_site() {
    // Positions (0-based): 1000 shared, 3000 new, 5000 weak in the parent, 7000 lost, 9000 weak in the child,
    // 11000 not covered in the parent.
    let parent = run_strain(&[site(1000, 30, 30), site(3000, 0, 30), site(5000, 4, 30), site(7000, 30, 30), site(9000, 30, 30)]);
    let child = run_strain(&[site(1000, 30, 30), site(3000, 30, 30), site(5000, 30, 30), site(7000, 0, 30), site(9000, 4, 30), site(11000, 30, 30)]);
    let inputs = [
        CombineInput { name: "parent", vcf: &parent.vcf, depth_bed: Some(&parent.bed), evidence_vcf: Some(&parent.evidence) },
        CombineInput { name: "child", vcf: &child.vcf, depth_bed: Some(&child.bed), evidence_vcf: Some(&child.evidence) },
    ];

    let c = combine(&inputs, &CombineOptions { parent: Some("parent"), min_depth: 5.0 }).unwrap();

    let status_at = |pos0: u64| c.rows.iter().find(|r| r.pos == pos0 + 1).map(|r| r.status[1]).unwrap_or_else(|| panic!("no row at {pos0}"));
    assert_eq!(status_at(1000), Status::Inherited);
    assert_eq!(status_at(3000), Status::New);
    assert_eq!(status_at(5000), Status::Unknown, "the parent has 4 of 30 alternate reads there");
    assert_eq!(status_at(7000), Status::Lost);
    assert_eq!(status_at(9000), Status::Unknown, "the child has 4 of 30 alternate reads there");
    assert_eq!(status_at(11000), Status::Unknown, "the parent has no reads there");
    // The weak evidence of the parent is shown as a weak cell, not as a plain `ref`.
    let weak = c.rows.iter().find(|r| r.pos == 5001).unwrap();
    assert!(matches!(weak.cells[0], Cell::Weak { depth: 30, .. }), "{:?}", weak.cells[0]);
}

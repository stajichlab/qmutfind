//! Confirming an introduced allele (a marker cassette that replaces or is inserted into the genome): where it joins the genome, whether
//! it is intact, and how many copies there are. The reference is the unmodified genome plus the cassette as an extra contig; the reads
//! come from the modified genome. Everything is simulated, so the right answer is known.

mod common;

use common::{reverse_complement, Lcg};
use qmutfind_core::{CheckStatus as Status, ConfirmOptions, Confirmer, Orientation, ReadKind, Target, Verdict};

const GENOME_LEN: usize = 120_000;
const CASSETTE_LEN: usize = 2_000;

struct World {
    genome: Vec<u8>,
    cassette: Vec<u8>,
}

fn world() -> World {
    let mut rng = Lcg(5);
    World { genome: (0..GENOME_LEN).map(|_| rng.base()).collect(), cassette: (0..CASSETTE_LEN).map(|_| rng.base()).collect() }
}

/// The genome with `replaced` (0-based, end exclusive) taken out and `insert` put in its place.
fn splice(genome: &[u8], replaced: std::ops::Range<usize>, insert: &[u8]) -> Vec<u8> {
    [&genome[..replaced.start], insert, &genome[replaced.end..]].concat()
}

/// FASTQ of 150-base reads from `seq`: uniform start positions, both strands, `depth` times coverage.
fn reads(seq: &[u8], depth: f64, seed: u64) -> String {
    let mut rng = Lcg(seed);
    let n = (seq.len() as f64 * depth / 150.0) as usize;
    let mut out = String::new();
    let mut state = seed.wrapping_mul(2862933555777941757).wrapping_add(3037000493);
    for i in 0..n {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let start = ((state >> 20) as usize) % (seq.len() - 150);
        let mut read = seq[start..start + 150].to_vec();
        if (state >> 40) & 1 == 1 {
            read = reverse_complement(&read);
        }
        let _ = rng.base();
        out.push_str(&format!("@r{i}\n{}\n+\n{}\n", String::from_utf8_lossy(&read), "I".repeat(150)));
    }
    out
}

fn reference_fasta(w: &World) -> Vec<u8> {
    format!(">chr\n{}\n>cass\n{}\n", String::from_utf8_lossy(&w.genome), String::from_utf8_lossy(&w.cassette)).into_bytes()
}

fn target() -> Target {
    Target { contig: "chr".into(), start: 50_001, end: 51_500 }
}

fn run(w: &World, fastq: &str, target: Option<Target>) -> qmutfind_core::Report {
    let options = ConfirmOptions { target, ..ConfirmOptions::default() };
    let mut c = Confirmer::new(&reference_fasta(w), &["cass"], ReadKind::Short, options).unwrap();
    for chunk in fastq.as_bytes().chunks(1 << 20) {
        c.add_reads(chunk).unwrap();
    }
    c.end_file().unwrap();
    c.report()
}

fn status(r: &qmutfind_core::Report, name: &str) -> Status {
    r.check(name).unwrap_or_else(|| panic!("no check {name}: {:?}", r.checks.iter().map(|c| c.name).collect::<Vec<_>>())).status
}

#[test]
fn an_intact_replacement_is_found_at_the_target_with_both_junctions_and_one_copy() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);

    let r = run(&w, &reads(&modified, 30.0, 1), Some(target()));

    assert_eq!(r.sites.len(), 1, "{:?}", r.sites);
    let site = &r.sites[0];
    assert_eq!(site.contig, "chr");
    assert_eq!((site.left, site.right), (Some(50_000), Some(51_501)), "the last genome base kept on the left and the first kept on the right");
    assert_eq!(site.orientation, Some(Orientation::Forward));
    assert_eq!(site.construct, "cass");
    assert!(site.reads_left >= 4 && site.reads_right >= 4, "{site:?}");
    assert_eq!(site.replaced_bases, Some(1_500));
    for name in ["integration_found", "both_sides_joined", "single_site", "at_target", "target_replaced", "cassette_coverage", "cassette_sequence", "cassette_structure", "copy_number"] {
        assert_eq!(status(&r, name), Status::Pass, "{name}: {}", r.check(name).unwrap().detail);
    }
    assert_eq!(r.verdict(), Verdict::Confirmed);
    assert!((r.cassettes[0].copy_ratio - 1.0).abs() < 0.3, "{}", r.cassettes[0].copy_ratio);
}

#[test]
fn a_cassette_that_went_in_the_other_way_is_reported_as_reverse() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &reverse_complement(&w.cassette));

    let r = run(&w, &reads(&modified, 30.0, 2), Some(target()));

    assert_eq!(r.sites.len(), 1, "{:?}", r.sites);
    assert_eq!((r.sites[0].left, r.sites[0].right), (Some(50_000), Some(51_501)));
    assert_eq!(r.sites[0].orientation, Some(Orientation::Reverse));
    assert_eq!(r.verdict(), Verdict::Confirmed, "{:?}", r.checks.iter().map(|c| (c.name, c.status)).collect::<Vec<_>>());
}

#[test]
fn a_point_change_and_a_deletion_inside_the_cassette_are_listed_and_fail_the_sequence_check() {
    let w = world();
    let mut damaged = w.cassette.clone();
    damaged[1_000] = if damaged[1_000] == b'A' { b'C' } else { b'A' }; // cassette position 1,001
    damaged.drain(1_500..1_530); // 30 bases gone after position 1,500
    let modified = splice(&w.genome, 50_000..51_500, &damaged);

    let r = run(&w, &reads(&modified, 30.0, 3), Some(target()));

    assert_eq!(status(&r, "cassette_sequence"), Status::Fail);
    let detail = &r.check("cassette_sequence").unwrap().detail;
    assert!(detail.contains("1001"), "{detail}");
    assert!(detail.contains("1500") || detail.contains("1501"), "{detail}");
    assert_eq!(status(&r, "integration_found"), Status::Pass);
    assert_eq!(status(&r, "both_sides_joined"), Status::Pass);
    assert_eq!(r.verdict(), Verdict::NotConfirmed);
    assert!(r.cassettes[0].variants.len() >= 2, "{:?}", r.cassettes[0].variants);
}

#[test]
fn a_missing_stretch_of_the_cassette_shows_as_a_coverage_gap_and_an_internal_junction() {
    let w = world();
    let mut broken = w.cassette.clone();
    broken.drain(800..1_200);
    let modified = splice(&w.genome, 50_000..51_500, &broken);

    let r = run(&w, &reads(&modified, 30.0, 4), Some(target()));

    assert_eq!(status(&r, "cassette_coverage"), Status::Fail);
    let detail = &r.check("cassette_coverage").unwrap().detail;
    assert!(detail.contains("801") || detail.contains("8"), "the gap is reported with its position: {detail}");
    assert_eq!(status(&r, "cassette_structure"), Status::Fail);
    assert!(r.check("cassette_structure").unwrap().detail.contains("cass"), "{}", r.check("cassette_structure").unwrap().detail);
    assert_eq!(r.verdict(), Verdict::NotConfirmed);
}

#[test]
fn a_second_copy_somewhere_else_is_found_as_a_second_site_and_raises_the_copy_number() {
    let w = world();
    let once = splice(&w.genome, 50_000..51_500, &w.cassette);
    let twice = splice(&once, 90_500..90_500, &w.cassette); // an extra insertion, nothing removed; 500 bases further on after the first insertion: reference position 90,000

    let r = run(&w, &reads(&twice, 30.0, 5), Some(target()));

    assert_eq!(r.sites.len(), 2, "{:?}", r.sites);
    assert_eq!(status(&r, "single_site"), Status::Fail);
    assert!(r.check("single_site").unwrap().detail.contains("90000"), "{}", r.check("single_site").unwrap().detail);
    assert_eq!(status(&r, "at_target"), Status::Pass, "the site at the target is still found");
    assert_eq!(status(&r, "copy_number"), Status::Fail);
    assert!((r.cassettes[0].copy_ratio - 2.0).abs() < 0.5, "{}", r.cassettes[0].copy_ratio);
    assert_eq!(r.verdict(), Verdict::NotConfirmed);
}

#[test]
fn a_population_that_is_half_wild_type_keeps_the_junctions_but_fails_the_replacement_and_copy_checks() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);
    let mixed = reads(&modified, 15.0, 6) + &reads(&w.genome, 15.0, 7).replace("@r", "@w");

    let r = run(&w, &mixed, Some(target()));

    assert_eq!(status(&r, "integration_found"), Status::Pass);
    assert_eq!(status(&r, "target_replaced"), Status::Fail, "{}", r.check("target_replaced").unwrap().detail);
    assert_eq!(status(&r, "copy_number"), Status::Fail, "{}", r.check("copy_number").unwrap().detail);
    assert_eq!(r.verdict(), Verdict::NotConfirmed);
}

#[test]
fn reads_of_the_unmodified_genome_find_no_integration_and_the_target_is_still_there() {
    let w = world();

    let r = run(&w, &reads(&w.genome, 30.0, 8), Some(target()));

    assert!(r.sites.is_empty());
    assert_eq!(status(&r, "integration_found"), Status::Fail);
    assert_eq!(status(&r, "target_replaced"), Status::Fail);
    assert_eq!(r.verdict(), Verdict::NotConfirmed);
}

#[test]
fn too_little_coverage_gives_inconclusive_not_a_pass_or_a_failure() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);

    let r = run(&w, &reads(&modified, 2.0, 9), Some(target()));

    assert_eq!(status(&r, "integration_found"), Status::Inconclusive, "{}", r.check("integration_found").unwrap().detail);
    assert_eq!(status(&r, "cassette_coverage"), Status::Inconclusive);
    assert_eq!(status(&r, "target_replaced"), Status::Inconclusive);
    assert_eq!(r.verdict(), Verdict::Inconclusive);
    assert!(r.check("integration_found").unwrap().detail.contains("coverage"), "{}", r.check("integration_found").unwrap().detail);
}

#[test]
fn an_insertion_away_from_the_target_is_not_at_the_target_and_leaves_the_target_in_place() {
    let w = world();
    let elsewhere = splice(&w.genome, 90_000..90_000, &w.cassette);

    let r = run(&w, &reads(&elsewhere, 30.0, 10), Some(target()));

    assert_eq!(r.sites.len(), 1);
    assert_eq!((r.sites[0].left, r.sites[0].right), (Some(90_000), Some(90_001)));
    assert_eq!(r.sites[0].replaced_bases, Some(0));
    assert_eq!(status(&r, "at_target"), Status::Fail);
    assert!(r.check("at_target").unwrap().detail.contains("90000"), "{}", r.check("at_target").unwrap().detail);
    assert_eq!(status(&r, "target_replaced"), Status::Fail);
    assert_eq!(r.verdict(), Verdict::NotConfirmed);
}

#[test]
fn without_a_target_the_target_checks_do_not_apply_and_the_site_is_still_reported() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);

    let r = run(&w, &reads(&modified, 30.0, 11), None);

    assert_eq!(status(&r, "at_target"), Status::NotApplicable);
    assert_eq!(status(&r, "target_replaced"), Status::NotApplicable);
    assert_eq!(r.sites.len(), 1);
    assert_eq!(r.verdict(), Verdict::Confirmed);
}

#[test]
fn the_text_report_names_every_check_and_the_junction_coordinates() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);
    let r = run(&w, &reads(&modified, 30.0, 12), Some(target()));

    let text = r.to_text();

    for name in ["integration_found", "at_target", "cassette_coverage", "copy_number"] {
        assert!(text.contains(name), "{text}");
    }
    assert!(text.contains("50000") && text.contains("51501"), "{text}");
    assert!(text.contains("CONFIRMED"), "{text}");
    assert!(text.contains("Junctions that involve the cassette"), "{text}");
    let json = r.to_json();
    assert!(json.contains("\"verdict\":\"confirmed\"") && json.contains("\"integration_found\""), "{json}");
}

/// A sweep over many random genomes and cassettes, both orientations, with the replaced length and the place varied. Slow (about a minute):
/// run it with `cargo test -p qmutfind-core --test confirm -- --ignored`.
#[test]
#[ignore]
fn intact_replacements_are_confirmed_with_the_exact_flanks_in_many_random_worlds() {
    let mut failures = Vec::new();
    for seed in 0..24u64 {
        let mut rng = Lcg(1000 + seed);
        let genome: Vec<u8> = (0..GENOME_LEN).map(|_| rng.base()).collect();
        let cassette: Vec<u8> = (0..CASSETTE_LEN).map(|_| rng.base()).collect();
        let w = World { genome, cassette };
        let start = 20_000 + (seed as usize * 3_137) % 80_000;
        let length = 300 + (seed as usize * 211) % 2_500;
        let reverse = seed % 2 == 1;
        let insert = if reverse { reverse_complement(&w.cassette) } else { w.cassette.clone() };
        let modified = splice(&w.genome, start..start + length, &insert);
        let t = Target { contig: "chr".into(), start: start as u64 + 1, end: (start + length) as u64 };
        let r = run(&w, &reads(&modified, 30.0, 500 + seed), Some(t));
        let ok = r.verdict() == Verdict::Confirmed
            && r.sites.len() == 1
            && r.sites[0].left == Some(start as u64)
            && r.sites[0].right == Some((start + length) as u64 + 1)
            && r.sites[0].orientation == Some(if reverse { Orientation::Reverse } else { Orientation::Forward });
        if !ok {
            failures.push(format!("seed {seed} (replace {length} at {start}, reverse {reverse}): verdict {:?}, sites {:?}, checks {:?}", r.verdict(), r.sites, r.checks.iter().filter(|c| c.status != Status::Pass).map(|c| (c.name, c.detail.clone())).collect::<Vec<_>>()));
        }
    }
    assert!(failures.is_empty(), "{} of 24 worlds failed:\n{}", failures.len(), failures.join("\n"));
}

#[test]
fn a_gene_is_found_by_its_id_name_or_locus_tag_in_a_refseq_style_gff() {
    let gff = "##gff-version 3\n\
        chr\tRefSeq\tgene\t216\t871\t.\t+\t.\tID=gene-AFUA_1G00100;Dbxref=GeneID:1;Name=AFUA_1G00100;gbkey=Gene;locus_tag=AFUA_1G00100;old_locus_tag=Afu1g00100\n\
        chr\tRefSeq\tmRNA\t216\t871\t.\t+\t.\tID=rna-XM_1;Parent=gene-AFUA_1G00100\n\
        chr\tRefSeq\tgene\t3752\t4495\t.\t-\t.\tID=gene-crpA;Name=crpA;gene=crpA;locus_tag=AFUA_1G00110\n";
    for name in ["AFUA_1G00100", "gene-AFUA_1G00100", "Afu1g00100"] {
        assert_eq!(qmutfind_core::find_gene(gff.as_bytes(), name).unwrap(), Target { contig: "chr".into(), start: 216, end: 871 }, "{name}");
    }
    assert_eq!(qmutfind_core::find_gene(gff.as_bytes(), "crpA").unwrap(), Target { contig: "chr".into(), start: 3752, end: 4495 });
    assert_eq!(qmutfind_core::find_gene(gff.as_bytes(), "AFUA_1G00110").unwrap().start, 3752);
    let err = qmutfind_core::find_gene(gff.as_bytes(), "XM_1").unwrap_err();
    assert!(err.contains("XM_1"), "an mRNA line is not a gene: {err}");
}

fn run_with_construct(w: &World, construct: &[u8], fastq: &str) -> qmutfind_core::Report {
    let fasta = format!(">chr\n{}\n>cass\n{}\n", String::from_utf8_lossy(&w.genome), String::from_utf8_lossy(construct)).into_bytes();
    let mut c = Confirmer::new(&fasta, &["cass"], ReadKind::Short, ConfirmOptions { target: Some(target()), ..ConfirmOptions::default() }).unwrap();
    c.add_reads(fastq.as_bytes()).unwrap();
    c.end_file().unwrap();
    c.report()
}

#[test]
fn a_construct_with_homology_arms_gets_a_warning_that_names_the_shared_stretches() {
    let w = world();
    // A full donor: 300 bases of the left arm, the cassette, 300 bases of the right arm (the arms are genome sequence).
    let donor = [&w.genome[49_700..50_000], &w.cassette[..], &w.genome[51_500..51_800]].concat();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);

    let r = run_with_construct(&w, &donor, &reads(&modified, 30.0, 13));

    assert!(!r.warnings.is_empty(), "no warning for a construct that holds genome sequence");
    let text = r.warnings.join(" ");
    assert!(text.contains("cass") && text.contains("1-300") || text.contains("1-"), "{text}");
    assert!(text.contains("arm"), "{text}");
    assert!(r.to_text().contains("Warning"), "{}", r.to_text());
    assert!(r.to_json().contains("\"warnings\""), "{}", r.to_json());
}

#[test]
fn a_cassette_without_arms_gets_no_warning() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);

    let r = run(&w, &reads(&modified, 30.0, 14), Some(target()));

    assert!(r.warnings.is_empty(), "{:?}", r.warnings);
}

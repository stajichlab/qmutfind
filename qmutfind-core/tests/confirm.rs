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
    // What the page draws from: the cassette (depth summary, gaps, variants), the junctions on the cassette, and the target.
    assert!(json.contains("\"cassettes\":[{\"name\":\"cass\",\"length\":2000"), "{json}");
    assert!(json.contains("\"gaps\":[]") && json.contains("\"variants\":[]"), "{json}");
    assert!(json.contains("\"construct_at_left\":1,") && json.contains("\"construct_at_right\":2000"), "{json}");
    assert!(json.contains("\"target\":{\"contig\":\"chr\",\"start\":50001,\"end\":51500}"), "{json}");
    assert!(json.contains("\"homology_left\":") && json.contains("\"copy_ratio\":"), "{json}");
}

#[test]
fn the_json_lists_gaps_and_variants_of_a_damaged_cassette_and_escapes_text() {
    let w = world();
    let mut damaged = w.cassette.clone();
    damaged[1_000] = if damaged[1_000] == b'A' { b'C' } else { b'A' };
    damaged.drain(800..1_200);
    let modified = splice(&w.genome, 50_000..51_500, &damaged);
    let r = run(&w, &reads(&modified, 30.0, 16), Some(target()));
    let json = r.to_json();
    assert!(json.contains("\"gaps\":[[8"), "the gap starts around 801: {json}");
    assert!(json.contains("\"verdict\":\"not_confirmed\""), "{json}");
    // Balanced braces and brackets outside strings: a cheap check that the JSON is well formed.
    let (mut depth, mut in_string, mut escaped) = (0i32, false, false);
    for ch in json.chars() {
        if in_string { if escaped { escaped = false } else if ch == '\\' { escaped = true } else if ch == '"' { in_string = false } continue }
        match ch { '"' => in_string = true, '{' | '[' => depth += 1, '}' | ']' => depth -= 1, _ => {} }
        assert!(depth >= 0, "{json}");
    }
    assert_eq!(depth, 0, "{json}");
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

#[test]
fn a_confirmer_can_be_made_from_a_genome_fasta_and_a_construct_fasta_and_checks_the_names() {
    let w = world();
    let genome = format!(">chr\n{}\n", String::from_utf8_lossy(&w.genome)).into_bytes();
    let construct = format!(">hph\n{}\n", String::from_utf8_lossy(&w.cassette)).into_bytes();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);
    let options = ConfirmOptions { target: Some(target()), ..ConfirmOptions::default() };

    let mut c = Confirmer::from_fasta(&genome, &construct, ReadKind::Short, options.clone()).unwrap();
    c.add_reads(reads(&modified, 30.0, 15).as_bytes()).unwrap();
    c.end_file().unwrap();
    let r = c.report();

    assert_eq!(r.verdict(), Verdict::Confirmed);
    assert_eq!(r.cassettes[0].name, "hph");
    // The tables that the CLI and the page write.
    let junctions = r.junctions_tsv();
    assert!(junctions.starts_with("contig_a\tpos_a\tside_a\tcontig_b\tpos_b\tside_b\treads\toverlap\tkind\n"), "{junctions}");
    assert!(junctions.lines().filter(|l| l.ends_with("\tintegration")).count() == 2, "{junctions}");
    let depth = r.depth_tsv();
    assert_eq!(depth.lines().count(), 1 + CASSETTE_LEN);
    assert!(depth.starts_with("contig\tpos\tdepth\nhph\t1\t"), "{}", &depth[..40]);

    let clash = format!(">chr\n{}\n", String::from_utf8_lossy(&w.cassette)).into_bytes();
    let err = Confirmer::from_fasta(&genome, &clash, ReadKind::Short, options.clone()).err().unwrap();
    assert!(err.contains("chr") && err.contains("same name"), "{err}");
    let err = Confirmer::from_fasta(&genome, b"", ReadKind::Short, options).err().unwrap();
    assert!(err.contains("construct") && err.contains("no sequence"), "{err}");
}

// ---- A full donor: homology arms around the cassette. The arms are cut off before mapping, and where they lie in the genome gives the target.

fn donor(w: &World) -> Vec<u8> {
    [&w.genome[49_700..50_000], &w.cassette[..], &w.genome[51_500..51_800]].concat()
}

fn confirm_donor(w: &World, donor: &[u8], fastq: &str, target: Option<Target>, strip_arms: bool) -> qmutfind_core::Report {
    let genome = format!(">chr\n{}\n", String::from_utf8_lossy(&w.genome)).into_bytes();
    let construct = format!(">donor\n{}\n", String::from_utf8_lossy(donor)).into_bytes();
    let mut c = Confirmer::from_fasta(&genome, &construct, ReadKind::Short, ConfirmOptions { target, strip_arms, ..ConfirmOptions::default() }).unwrap();
    c.add_reads(fastq.as_bytes()).unwrap();
    c.end_file().unwrap();
    c.report()
}

#[test]
fn homology_arms_are_cut_off_the_target_is_inferred_and_the_cassette_is_confirmed() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);

    let r = confirm_donor(&w, &donor(&w), &reads(&modified, 30.0, 20), None, true);

    assert_eq!((r.cassettes[0].trimmed_left, r.cassettes[0].trimmed_right), (300, 300));
    assert_eq!(r.cassettes[0].length, 2_000, "the mapped construct is the cassette without its arms");
    assert!(r.target_inferred);
    assert_eq!(r.target, Some(Target { contig: "chr".into(), start: 50_001, end: 51_500 }));
    assert_eq!((r.sites[0].left, r.sites[0].right), (Some(50_000), Some(51_501)));
    assert_eq!(r.sites[0].construct_at_left, Some(1), "positions are on the cassette without the arm");
    assert_eq!(r.verdict(), Verdict::Confirmed, "{:?}", r.checks.iter().map(|c| (c.name, c.status, c.detail.clone())).collect::<Vec<_>>());
    assert!(r.warnings.is_empty(), "{:?}", r.warnings);
    assert!(r.notes.iter().any(|n| n.contains("arm") && n.contains("300")), "{:?}", r.notes);
    assert!(r.to_text().contains("Note:") && r.to_text().contains("inferred"), "{}", r.to_text());
}

#[test]
fn a_donor_in_the_other_orientation_gives_the_same_target_and_a_reverse_cassette() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);

    let r = confirm_donor(&w, &reverse_complement(&donor(&w)), &reads(&modified, 30.0, 21), None, true);

    assert_eq!(r.target, Some(Target { contig: "chr".into(), start: 50_001, end: 51_500 }), "{:?}", r.notes);
    assert_eq!(r.sites[0].orientation, Some(Orientation::Reverse), "the cassette lies the other way round compared with the construct file");
    assert_eq!(r.verdict(), Verdict::Confirmed);
}

#[test]
fn a_target_that_agrees_with_the_arms_is_kept_and_not_called_inferred() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);

    let r = confirm_donor(&w, &donor(&w), &reads(&modified, 30.0, 22), Some(target()), true);

    assert!(!r.target_inferred);
    assert!(r.warnings.is_empty(), "{:?}", r.warnings);
    assert_eq!(r.verdict(), Verdict::Confirmed);
}

#[test]
fn a_target_that_the_arms_contradict_gets_a_warning_and_fails_at_target() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);
    let wrong = Target { contig: "chr".into(), start: 80_001, end: 81_500 };

    let r = confirm_donor(&w, &donor(&w), &reads(&modified, 30.0, 23), Some(wrong), true);

    assert!(r.warnings.iter().any(|w| w.contains("arms") && w.contains("50001") && w.contains("51500")), "{:?}", r.warnings);
    assert_eq!(status(&r, "at_target"), Status::Fail);
    assert_eq!(r.verdict(), Verdict::NotConfirmed);
}

#[test]
fn a_donor_with_one_arm_is_trimmed_on_that_side_and_no_target_is_inferred() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);
    let left_only = [&w.genome[49_700..50_000], &w.cassette[..]].concat();

    let r = confirm_donor(&w, &left_only, &reads(&modified, 30.0, 24), None, true);

    assert_eq!((r.cassettes[0].trimmed_left, r.cassettes[0].trimmed_right), (300, 0));
    assert!(!r.target_inferred && r.target.is_none());
    assert!(r.notes.iter().any(|n| n.contains("one arm") || n.contains("not inferred")), "{:?}", r.notes);
    assert_eq!(r.sites.len(), 1);
    assert_eq!(status(&r, "at_target"), Status::NotApplicable);
}

#[test]
fn arms_can_be_kept_and_then_the_old_warning_stays() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);

    let r = confirm_donor(&w, &donor(&w), &reads(&modified, 30.0, 25), None, false);

    assert_eq!(r.cassettes[0].length, 2_600);
    assert!(r.warnings.iter().any(|w| w.contains("arm")), "{:?}", r.warnings);
}

#[test]
fn genome_sequence_inside_the_construct_is_not_cut_but_warned_about() {
    let w = world();
    // 400 genome bases in the middle of the cassette (not at an end): they cannot be arms.
    let inner = [&w.cassette[..1_000], &w.genome[10_000..10_400], &w.cassette[1_000..]].concat();
    let modified = splice(&w.genome, 50_000..51_500, &inner);

    let r = confirm_donor(&w, &inner, &reads(&modified, 30.0, 26), None, true);

    assert_eq!((r.cassettes[0].trimmed_left, r.cassettes[0].trimmed_right), (0, 0));
    assert!(r.warnings.iter().any(|w| w.contains("also occur in the genome")), "{:?}", r.warnings);
}

/// Donors with arms of 60 to 800 bases in both orientations, at many places: the arms are cut off (to within a few bases that the cassette shares with the
/// genome by chance), the target is inferred to within the same bases, and the verdict is confirmed. Slow: run it with `-- --ignored`.
#[test]
#[ignore]
fn donors_with_arms_are_trimmed_and_the_target_is_inferred_in_many_random_worlds() {
    let mut failures = Vec::new();
    for seed in 0..16u64 {
        let mut rng = Lcg(2000 + seed);
        let genome: Vec<u8> = (0..GENOME_LEN).map(|_| rng.base()).collect();
        let cassette: Vec<u8> = (0..CASSETTE_LEN).map(|_| rng.base()).collect();
        let w = World { genome, cassette };
        let start = 20_000 + (seed as usize * 4_111) % 70_000;
        let length = 300 + (seed as usize * 313) % 2_000;
        let (arm_left, arm_right) = (60 + (seed as usize * 47) % 740, 60 + (seed as usize * 91) % 740);
        let donor = [&w.genome[start - arm_left..start], &w.cassette[..], &w.genome[start + length..start + length + arm_right]].concat();
        let given = if seed % 2 == 1 { reverse_complement(&donor) } else { donor };
        let modified = splice(&w.genome, start..start + length, &w.cassette);
        let r = confirm_donor(&w, &given, &reads(&modified, 30.0, 700 + seed), None, true);
        let (want_left, want_right) = if seed % 2 == 1 { (arm_right, arm_left) } else { (arm_left, arm_right) };
        // The edge of an arm is known up to the bases that the cassette shares with the genome by chance (about a quarter of the time one base, an eighth of
        // the time two): the arm may be cut that many bases too long, and the target that many bases too short at each end.
        let slack = 6u64;
        let (tl, tr) = (r.cassettes[0].trimmed_left, r.cassettes[0].trimmed_right);
        let target_ok = match &r.target {
            Some(t) => t.contig == "chr" && t.start >= start as u64 + 1 && t.start <= start as u64 + 1 + slack && t.end <= (start + length) as u64 && t.end + slack >= (start + length) as u64,
            None => false,
        };
        let ok = tl >= want_left as u64 && tl <= want_left as u64 + slack && tr >= want_right as u64 && tr <= want_right as u64 + slack && target_ok && r.verdict() == Verdict::Confirmed;
        if !ok {
            failures.push(format!("seed {seed} (arms {arm_left}/{arm_right}, replace {length} at {start}, reversed {}): trimmed {:?}, target {:?}, verdict {:?}, notes {:?}", seed % 2 == 1, (r.cassettes[0].trimmed_left, r.cassettes[0].trimmed_right), r.target, r.verdict(), r.notes));
        }
    }
    assert!(failures.is_empty(), "{} of 16 worlds failed:\n{}", failures.len(), failures.join("\n"));
}

// ---- Long reads (Nanopore-like): 3 to 9 kb, about 8% errors (substitutions and insertions and deletions), both strands.

fn long_reads(seq: &[u8], depth: f64, seed: u64) -> String {
    let mut state = seed.wrapping_mul(2862933555777941757).wrapping_add(3037000493);
    let mut next = move || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 20) as usize
    };
    let mut out = String::new();
    let mut total = 0usize;
    let mut i = 0;
    while (total as f64) < seq.len() as f64 * depth {
        let len = 3_000 + next() % 6_000;
        let start = next() % (seq.len() - len);
        let mut read = Vec::with_capacity(len + 200);
        for &b in &seq[start..start + len] {
            match next() % 1000 {
                0..=29 => read.push(match b { b'A' => b'C', b'C' => b'G', b'G' => b'T', _ => b'A' }), // substitution 3%
                30..=54 => { read.push(b); read.push(b"ACGT"[next() % 4]); }                          // insertion 2.5%
                55..=79 => {}                                                                          // deletion 2.5%
                _ => read.push(b),
            }
        }
        if next() & 1 == 1 {
            read = reverse_complement(&read);
        }
        total += len;
        out.push_str(&format!("@l{i}\n{}\n+\n{}\n", String::from_utf8_lossy(&read), "5".repeat(read.len())));
        i += 1;
    }
    out
}

fn run_long(w: &World, fastq: &str, target: Option<Target>) -> qmutfind_core::Report {
    let options = ConfirmOptions { target, ..ConfirmOptions::defaults_for(ReadKind::Long) };
    let mut c = Confirmer::new(&reference_fasta(w), &["cass"], ReadKind::Long, options).unwrap();
    for chunk in fastq.as_bytes().chunks(1 << 20) {
        c.add_reads(chunk).unwrap();
    }
    c.end_file().unwrap();
    c.report()
}

fn within(value: Option<u64>, truth: u64, slack: u64) -> bool {
    value.map_or(false, |v| v.abs_diff(truth) <= slack)
}

#[test]
#[cfg_attr(debug_assertions, ignore = "slow in a debug build: run with --release (CI does)")]
fn long_reads_confirm_an_intact_replacement_and_show_that_whole_reads_span_the_insertion() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);

    let r = run_long(&w, &long_reads(&modified, 15.0, 40), Some(target()));

    assert_eq!(r.sites.len(), 1, "{:?}", r.sites);
    assert!(within(r.sites[0].left, 50_000, 20) && within(r.sites[0].right, 51_501, 20), "long-read breakpoints are uncertain by a few bases: {:?}", r.sites);
    assert_eq!(r.sites[0].orientation, Some(Orientation::Forward));
    assert!(r.spanning.reads >= 3 && r.spanning.full >= 3, "{:?}", r.spanning);
    assert_eq!(status(&r, "whole_insertion_in_one_read"), Status::Pass, "{}", r.check("whole_insertion_in_one_read").unwrap().detail);
    assert_eq!(r.verdict(), Verdict::Confirmed, "{:?}", r.checks.iter().filter(|c| c.status != Status::Pass).map(|c| (c.name, c.status, c.detail.clone())).collect::<Vec<_>>());
}

#[test]
#[cfg_attr(debug_assertions, ignore = "slow in a debug build: run with --release (CI does)")]
fn long_reads_show_a_missing_stretch_of_the_cassette_and_no_read_spans_it_whole() {
    let w = world();
    let mut broken = w.cassette.clone();
    broken.drain(800..1_200);
    let modified = splice(&w.genome, 50_000..51_500, &broken);

    let r = run_long(&w, &long_reads(&modified, 15.0, 41), Some(target()));

    assert_eq!(r.verdict(), Verdict::NotConfirmed, "{:?}", r.checks.iter().map(|c| (c.name, c.status)).collect::<Vec<_>>());
    assert!(status(&r, "cassette_coverage") == Status::Fail || status(&r, "cassette_structure") == Status::Fail || status(&r, "cassette_sequence") == Status::Fail);
    assert_ne!(status(&r, "whole_insertion_in_one_read"), Status::Pass, "reads span the insertion but not the whole cassette");
}

#[test]
#[cfg_attr(debug_assertions, ignore = "slow in a debug build: run with --release (CI does)")]
fn long_reads_find_the_wild_type_and_a_second_copy_like_short_reads_do() {
    let w = world();
    let r = run_long(&w, &long_reads(&w.genome, 15.0, 42), Some(target()));
    assert_eq!(status(&r, "integration_found"), Status::Fail);
    assert_eq!(status(&r, "target_replaced"), Status::Fail);

    let once = splice(&w.genome, 50_000..51_500, &w.cassette);
    let twice = splice(&once, 90_500..90_500, &w.cassette);
    let r = run_long(&w, &long_reads(&twice, 15.0, 43), Some(target()));
    assert_eq!(r.sites.len(), 2, "{:?}", r.sites);
    assert_eq!(status(&r, "single_site"), Status::Fail);
}

#[test]
fn for_short_reads_the_whole_insertion_check_does_not_apply() {
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);

    let r = run(&w, &reads(&modified, 30.0, 44), Some(target()));

    assert_eq!(status(&r, "whole_insertion_in_one_read"), Status::NotApplicable);
    assert!(r.check("whole_insertion_in_one_read").unwrap().detail.contains("short"), "{}", r.check("whole_insertion_in_one_read").unwrap().detail);
}

/// Long reads in many random worlds, both orientations, 15x: no world may be NotConfirmed or have a wrong site, and at most 20% may be Inconclusive (the depth
/// of long reads is noisy at this coverage). Slow in a debug build: `cargo test --release -p qmutfind-core --test confirm -- --ignored long_reads_in_many`.
#[test]
#[ignore]
fn long_reads_in_many_random_worlds_confirm_an_intact_replacement() {
    let mut failures = Vec::new();
    let mut inconclusive = 0;
    for seed in 0..40u64 {
        let mut rng = Lcg(3000 + seed);
        let genome: Vec<u8> = (0..GENOME_LEN).map(|_| rng.base()).collect();
        let cassette: Vec<u8> = (0..CASSETTE_LEN).map(|_| rng.base()).collect();
        let w = World { genome, cassette };
        let start = 20_000 + (seed as usize * 5_003) % 70_000;
        let length = 300 + (seed as usize * 277) % 2_500;
        let reverse = seed % 2 == 1;
        let insert = if reverse { reverse_complement(&w.cassette) } else { w.cassette.clone() };
        let modified = splice(&w.genome, start..start + length, &insert);
        let t = Target { contig: "chr".into(), start: start as u64 + 1, end: (start + length) as u64 };
        let r = run_long(&w, &long_reads(&modified, 15.0, 900 + seed), Some(t));
        let site_ok = r.sites.len() == 1
            && within(r.sites[0].left, start as u64, 20)
            && within(r.sites[0].right, (start + length) as u64 + 1, 20)
            && r.sites[0].orientation == Some(if reverse { Orientation::Reverse } else { Orientation::Forward });
        if r.verdict() == Verdict::Inconclusive {
            inconclusive += 1;
        }
        if !(site_ok && r.verdict() != Verdict::NotConfirmed) {
            failures.push(format!("seed {seed} (replace {length} at {start}, reverse {reverse}): verdict {:?}, sites {:?}, spanning {:?}, not passing {:?}", r.verdict(), r.sites, r.spanning, r.checks.iter().filter(|c| c.status != Status::Pass).map(|c| (c.name, c.status, c.detail.clone())).collect::<Vec<_>>()));
        }
    }
    println!("inconclusive: {inconclusive} of 40");
    assert!(failures.is_empty(), "{} of 40 worlds failed:\n{}", failures.len(), failures.join("\n"));
    assert!(inconclusive <= 8, "{inconclusive} of 40 worlds were inconclusive");
}

#[test]
#[cfg_attr(debug_assertions, ignore = "slow in a debug build: run with --release (CI does)")]
fn long_reads_confirm_when_the_cassette_is_much_longer_than_the_replaced_region() {
    // 300 genome bases replaced by 2,000 cassette bases: the aligner puts a read across this into one alignment with a long insertion.
    let w = world();
    let modified = splice(&w.genome, 20_000..20_300, &w.cassette);
    let t = Target { contig: "chr".into(), start: 20_001, end: 20_300 };

    let r = run_long(&w, &long_reads(&modified, 15.0, 50), Some(t));

    assert_eq!(r.sites.len(), 1, "{:?}", r.sites);
    assert!(within(r.sites[0].left, 20_000, 20) && within(r.sites[0].right, 20_301, 20), "{:?}", r.sites);
    assert!(r.spanning.full >= 2, "{:?}", r.spanning);
    assert_eq!(r.verdict(), Verdict::Confirmed, "{:?}", r.checks.iter().filter(|c| c.status != Status::Pass).map(|c| (c.name, c.status, c.detail.clone())).collect::<Vec<_>>());
}

#[test]
#[cfg_attr(debug_assertions, ignore = "slow in a debug build: run with --release (CI does)")]
fn long_reads_confirm_when_the_replaced_region_is_much_longer_than_the_cassette() {
    // 4,500 genome bases replaced by the 2,000 base cassette: the alignment has a long deletion.
    let w = world();
    let modified = splice(&w.genome, 30_000..34_500, &w.cassette);
    let t = Target { contig: "chr".into(), start: 30_001, end: 34_500 };

    let r = run_long(&w, &long_reads(&modified, 15.0, 51), Some(t));

    assert_eq!(r.sites.len(), 1, "{:?}", r.sites);
    assert!(within(r.sites[0].left, 30_000, 20) && within(r.sites[0].right, 34_501, 20), "{:?}", r.sites);
    assert!(r.spanning.full >= 2, "{:?}", r.spanning);
    assert_eq!(r.verdict(), Verdict::Confirmed, "{:?}", r.checks.iter().filter(|c| c.status != Status::Pass).map(|c| (c.name, c.status, c.detail.clone())).collect::<Vec<_>>());
}

#[test]
#[cfg_attr(debug_assertions, ignore = "slow in a debug build: run with --release (CI does)")]
fn long_reads_show_a_deletion_without_a_cassette_as_a_junction_between_genome_positions() {
    // The gene is deleted and nothing is inserted: the long reads cross the deletion in one alignment with a long deletion.
    let w = world();
    let deleted = splice(&w.genome, 40_000..42_500, &[]);
    let t = Target { contig: "chr".into(), start: 40_001, end: 42_500 };

    let r = run_long(&w, &long_reads(&deleted, 15.0, 52), Some(t));

    assert_eq!(status(&r, "target_replaced"), Status::Pass, "{}", r.check("target_replaced").unwrap().detail);
    assert_eq!(status(&r, "integration_found"), Status::Fail);
    assert!(r.junctions.iter().any(|j| j.a.contig == "chr" && j.b.contig == "chr" && within(Some(j.a.pos), 40_000, 20) && within(Some(j.b.pos), 42_501, 20)), "{:?}", r.junctions);
}

#[test]
#[cfg_attr(debug_assertions, ignore = "slow in a debug build: run with --release (CI does)")]
fn long_reads_with_half_wild_type_reads_give_an_unclear_copy_number_not_a_pass() {
    // At long-read coverage a half-and-half mixture (cassette depth about half of the genome depth) cannot be told from depth noise: it is unclear, never a pass.
    let w = world();
    let modified = splice(&w.genome, 50_000..51_500, &w.cassette);
    let mixed = long_reads(&modified, 8.0, 53) + &long_reads(&w.genome, 8.0, 54).replace("@l", "@w");

    let r = run_long(&w, &mixed, Some(target()));

    assert_ne!(status(&r, "copy_number"), Status::Pass, "{}", r.check("copy_number").unwrap().detail);
    assert_eq!(status(&r, "target_replaced"), Status::Fail, "{}", r.check("target_replaced").unwrap().detail);
    assert_eq!(r.verdict(), Verdict::NotConfirmed);
}

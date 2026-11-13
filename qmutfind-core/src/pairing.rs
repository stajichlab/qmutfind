//! Pair-aware choice of the alignments of the two mates of a fragment.
//!
//! Each mate is aligned on its own and gets a list of candidates: its primary alignments and its secondary ones. A mate in a
//! repeat has several equally good candidates and a mapping quality of 0, so on its own it is ignored by the caller. Its mate,
//! when it lies in unique sequence, says which candidate is right: the one that makes a proper pair with it. This is the rule
//! (the same idea as in `bwa mem` and `minimap2 -ax sr`, in a simpler form):
//!
//! - A proper pair is two alignments on one contig on opposite strands, the forward one starting first (it may not start
//!   after the reverse one), with at most `max_insert` bases from the start of the forward mate to the end of the reverse mate.
//! - If there is a proper pair among the candidates, the pair with the best total score is chosen (ties: both primary, then the
//!   shorter insert).
//! - If that pair is the only proper pair, each mate gets at least the mapping quality of the other (up to 60): a read that is
//!   ambiguous on its own is as reliable as the mate that fixes its place. If several pairs are proper, the choice is made but the
//!   qualities stay, because the repeat is not resolved.
//! - If there is no proper pair, nothing changes: every mate keeps its primary alignments.

use crate::Alignment;

/// The default largest insert: from the start of the forward mate to the end of the reverse mate, in bases.
pub const DEFAULT_MAX_INSERT: usize = 1000;

/// Two alignments make a proper pair: one contig, opposite strands, forward mate first, within `max_insert` bases.
fn is_proper(a: &Alignment, b: &Alignment, max_insert: usize) -> bool {
    if a.contig != b.contig || a.reverse == b.reverse {
        return false;
    }
    let (forward, reverse) = if a.reverse { (b, a) } else { (a, b) };
    forward.ref_start <= reverse.ref_start && reverse.ref_end() > forward.ref_start && reverse.ref_end() - forward.ref_start <= max_insert
}

fn insert_size(a: &Alignment, b: &Alignment) -> usize {
    a.ref_start.min(b.ref_start).abs_diff(a.ref_end().max(b.ref_end()))
}

/// The alignments to count for each mate. `candidates1` and `candidates2` are the candidates of the first and second mate
/// (`Aligner::align_candidates`). Returns the alignments of the first and of the second mate.
pub fn choose_pair(candidates1: Vec<Alignment>, candidates2: Vec<Alignment>, max_insert: usize) -> (Vec<Alignment>, Vec<Alignment>) {
    let primaries = |c: &[Alignment]| c.iter().filter(|a| a.is_primary()).cloned().collect::<Vec<_>>();
    let pairs: Vec<(usize, usize)> = candidates1
        .iter()
        .enumerate()
        .flat_map(|(i, a)| candidates2.iter().enumerate().filter(move |(_, b)| is_proper(a, b, max_insert)).map(move |(j, _)| (i, j)))
        .collect();
    let Some(&(i, j)) = pairs.iter().max_by_key(|&&(i, j)| {
        let (a, b) = (&candidates1[i], &candidates2[j]);
        (a.score() + b.score(), usize::from(a.is_primary()) + usize::from(b.is_primary()), std::cmp::Reverse(insert_size(a, b)), std::cmp::Reverse(i + j))
    }) else {
        return (primaries(&candidates1), primaries(&candidates2));
    };
    let (q1, q2) = (candidates1[i].mapq(), candidates2[j].mapq());
    let (n1, n2) = if pairs.len() == 1 { (q1.max(q2.min(60)), q2.max(q1.min(60))) } else { (q1, q2) };
    let keep = |candidates: Vec<Alignment>, chosen: usize, mapq: u8| -> Vec<Alignment> {
        if candidates[chosen].is_primary() {
            // The primary alignment stays, with the new quality; other primary alignments of a split read stay as they were.
            candidates
                .into_iter()
                .enumerate()
                .filter(|(k, a)| *k == chosen || a.is_primary())
                .map(|(k, a)| if k == chosen { a.with_mapq(mapq) } else { a })
                .collect()
        } else {
            // The mate says that the read belongs in a secondary place: that alignment replaces the primary one.
            vec![candidates.into_iter().nth(chosen).expect("chosen is an index of the candidates").with_mapq(mapq)]
        }
    };
    (keep(candidates1, i, n1), keep(candidates2, j, n2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::CigarOp;

    /// A 100-base alignment on `contig`.
    fn aln(contig: &str, start: usize, reverse: bool, mapq: u8, score: i32, primary: bool) -> Alignment {
        Alignment::new(contig, start, reverse, vec![CigarOp::Match(100)], 100, mapq, score, primary)
    }

    fn starts(v: &[Alignment]) -> Vec<(usize, u8)> {
        v.iter().map(|a| (a.ref_start, a.mapq())).collect()
    }

    #[test]
    fn a_mate_with_two_equal_candidates_takes_the_one_that_makes_a_proper_pair_and_the_quality_of_its_unique_mate() {
        let mate1 = vec![aln("c", 19_850, false, 60, 200, true)];
        let mate2 = vec![aln("c", 60_250, true, 0, 200, true), aln("c", 20_250, true, 0, 200, false)]; // the primary is the wrong copy

        let (a, b) = choose_pair(mate1, mate2, 1000);

        assert_eq!(starts(&a), vec![(19_850, 60)]);
        assert_eq!(starts(&b), vec![(20_250, 60)]); // the secondary alignment replaced the primary one and got quality 60
    }

    #[test]
    fn two_unique_mates_in_a_proper_pair_keep_their_own_qualities() {
        let (a, b) = choose_pair(vec![aln("c", 1000, false, 50, 190, true)], vec![aln("c", 1300, true, 40, 180, true)], 1000);

        assert_eq!((a[0].mapq(), b[0].mapq()), (50, 60.min(50).max(40))); // each mate gets at least the quality of the other
    }

    #[test]
    fn mates_that_are_not_a_proper_pair_keep_their_primary_alignments_and_qualities() {
        let m1 = || vec![aln("c", 1000, false, 60, 200, true)];
        let far = vec![aln("c", 9000, true, 0, 200, true), aln("c", 40_000, true, 0, 200, false)]; // 8,100 bases away
        let same_strand = vec![aln("c", 1300, false, 0, 200, true), aln("c", 5000, true, 0, 200, false)];
        let other_contig = vec![aln("d", 1300, true, 0, 200, true)];
        let dovetail = vec![aln("c", 900, true, 0, 200, true)]; // the reverse mate starts before the forward mate

        for (name, mate2) in [("far", far), ("same strand", same_strand), ("other contig", other_contig), ("dovetail", dovetail)] {
            let expected = starts(&mate2.iter().filter(|a| a.is_primary()).cloned().collect::<Vec<_>>());
            let (a, b) = choose_pair(m1(), mate2, 1000);
            assert_eq!(starts(&a), vec![(1000, 60)], "{name}");
            assert_eq!(starts(&b), expected, "{name}");
        }
    }

    #[test]
    fn when_both_copies_make_a_proper_pair_the_repeat_is_not_resolved_and_the_qualities_stay_at_zero() {
        // Both mates are in a two-copy repeat: (copy A, copy A) and (copy B, copy B) are both proper pairs.
        let mate1 = vec![aln("c", 1000, false, 0, 200, true), aln("c", 5000, false, 0, 200, false)];
        let mate2 = vec![aln("c", 1300, true, 0, 200, true), aln("c", 5300, true, 0, 200, false)];

        let (a, b) = choose_pair(mate1, mate2, 1000);

        assert_eq!((a[0].mapq(), b[0].mapq()), (0, 0));
        assert_eq!((a[0].ref_start, b[0].ref_start), (1000, 1300)); // the pair with both primaries is chosen
    }

    #[test]
    fn the_best_scoring_proper_pair_wins_and_with_two_proper_pairs_the_qualities_are_not_raised() {
        let mate1 = vec![aln("c", 1000, false, 60, 200, true)];
        let mate2 = vec![aln("c", 1300, true, 30, 150, true), aln("c", 1500, true, 0, 190, false)];

        let (a, b) = choose_pair(mate1, mate2, 1000);

        assert_eq!(starts(&b), vec![(1500, 0)]); // 200 + 190 beats 200 + 150; two proper pairs, so the quality of the chosen alignment stays
        assert_eq!(a[0].mapq(), 60);
    }

    #[test]
    fn a_mate_without_candidates_gives_no_pairing_and_the_other_mate_keeps_its_primary_alignments() {
        let (a, b) = choose_pair(vec![aln("c", 1000, false, 60, 200, true), aln("c", 8000, false, 0, 150, false)], Vec::new(), 1000);

        assert_eq!(starts(&a), vec![(1000, 60)]);
        assert!(b.is_empty());
    }

    #[test]
    fn the_insert_limit_is_a_setting() {
        let mate1 = || vec![aln("c", 1000, false, 60, 200, true)];
        let mate2 = || vec![aln("c", 1700, true, 0, 200, true)]; // insert 800

        assert_eq!(choose_pair(mate1(), mate2(), 1000).1[0].mapq(), 60);
        assert_eq!(choose_pair(mate1(), mate2(), 500).1[0].mapq(), 0);
    }
}

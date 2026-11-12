//! Variant size classes: SNP, small indel (1 to 49 bp) and large indel (50 bp and over).

use qmutfind_core::{classify, VariantClass};

const LARGE: usize = 50;

fn bases(n: usize) -> String {
    "A".repeat(n)
}

#[test]
fn a_single_base_change_is_a_snp() {
    assert_eq!(classify("A", "G", LARGE), VariantClass::Snp);
}

#[test]
fn an_indel_of_49_bp_is_small_and_of_50_bp_is_large_for_deletions_and_insertions() {
    // VCF form: one anchor base plus the changed bases.
    let del = |n: usize| classify(&format!("C{}", bases(n)), "C", LARGE);
    let ins = |n: usize| classify("C", &format!("C{}", bases(n)), LARGE);

    assert_eq!(del(1), VariantClass::SmallIndel);
    assert_eq!(del(49), VariantClass::SmallIndel);
    assert_eq!(del(50), VariantClass::LargeIndel);
    assert_eq!(ins(1), VariantClass::SmallIndel);
    assert_eq!(ins(49), VariantClass::SmallIndel);
    assert_eq!(ins(50), VariantClass::LargeIndel);
}

#[test]
fn the_large_indel_cutoff_is_a_setting() {
    assert_eq!(classify("CAAAAA", "C", 5), VariantClass::LargeIndel);
    assert_eq!(classify("CAAAAA", "C", 6), VariantClass::SmallIndel);
}

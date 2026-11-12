//! Left-alignment of indels, so calls from different tools agree on the position.

use qmutfind_core::{left_align_calls, Call};

fn call(pos: u64, reference: &str, alt: &str) -> Call {
    Call { contig: "chrA".to_string(), pos, reference: reference.to_string(), alt: alt.to_string(), depth: 20, alt_count: 19, svlen: None, qual: None, depth_ratio: None }
}

fn refs(seq: &str) -> Vec<(String, Vec<u8>)> {
    vec![("chrA".to_string(), seq.as_bytes().to_vec())]
}

fn shape(c: &Call) -> (u64, &str, &str) {
    (c.pos, c.reference.as_str(), c.alt.as_str())
}

#[test]
fn a_deletion_in_a_homopolymer_moves_to_the_base_before_the_run() {
    // CCCGTAAAACGT: the A run is positions 6 to 9. Deleting one A is written with the last two As.
    let out = left_align_calls(vec![call(8, "AA", "A")], &refs("CCCGTAAAACGT"));

    assert_eq!(shape(&out[0]), (5, "TA", "T"));
}

#[test]
fn an_insertion_after_a_dinucleotide_repeat_moves_to_the_base_before_the_repeat() {
    // TTGCACAGG: CA at 4-5 and 6-7. Inserting CA after position 7 equals inserting it after G at 3.
    let out = left_align_calls(vec![call(7, "A", "ACA")], &refs("TTGCACAGG"));

    assert_eq!(shape(&out[0]), (3, "G", "GCA"));
}

#[test]
fn snps_symbolic_calls_and_already_leftmost_indels_do_not_change() {
    let mut symbolic = call(7, "A", "<INS>");
    symbolic.svlen = Some(80);
    let input = vec![call(4, "C", "T"), symbolic.clone(), call(5, "TA", "T")];

    let out = left_align_calls(input.clone(), &refs("CCCGTAAAACGT"));

    assert_eq!(out, input);
}

#[test]
fn an_indel_that_cannot_be_anchored_to_the_left_is_left_as_it_was() {
    // Deleting the first A of AAAC at the very start of the contig has no base to its left.
    let input = vec![call(1, "AA", "A")];

    let out = left_align_calls(input.clone(), &refs("AAAC"));

    assert_eq!(out, input);
}

//! FASTQ parsing from byte chunks that can split a record anywhere, as a browser delivers them.

use qmutfind_core::{FastqChunker, FastqReader};

const FASTQ: &str = "@r1 extra\nACGT\n+\nIIII\n@r2\nGGGGCC\n+\nIIIIII\n@r3\nTTA\n+\nIII\n";

fn expected() -> Vec<(String, Vec<u8>)> {
    FastqReader::new(FASTQ.as_bytes()).map(|r| r.unwrap()).collect()
}

fn parse_in_chunks(text: &str, size: usize) -> std::io::Result<Vec<(String, Vec<u8>)>> {
    let mut chunker = FastqChunker::new();
    let mut out = Vec::new();
    for chunk in text.as_bytes().chunks(size) {
        out.extend(chunker.push(chunk)?);
    }
    out.extend(chunker.finish()?);
    Ok(out)
}

#[test]
fn the_records_are_the_same_for_every_chunk_size_including_one_byte_at_a_time() {
    for size in [1, 2, 3, 5, 7, 11, 13, 20, 1000] {
        assert_eq!(parse_in_chunks(FASTQ, size).unwrap(), expected(), "chunk size {size}");
    }
}

#[test]
fn a_last_record_without_a_trailing_newline_is_still_returned() {
    let text = FASTQ.trim_end_matches('\n');

    assert_eq!(parse_in_chunks(text, 4).unwrap(), expected());
}

#[test]
fn a_truncated_last_record_is_an_error_at_finish() {
    let truncated = "@r1\nACGT\n+\nIIII\n@r2\nGGGG\n";

    assert!(parse_in_chunks(truncated, 6).is_err());
}

#[test]
fn nothing_is_returned_until_a_whole_record_has_arrived() {
    let mut chunker = FastqChunker::new();

    assert!(chunker.push(b"@r1\nACGT\n+\nII").unwrap().is_empty());
    assert_eq!(chunker.push(b"II\n").unwrap(), vec![("r1".to_string(), b"ACGT".to_vec())]);
}

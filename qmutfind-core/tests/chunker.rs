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

#[test]
fn a_record_that_never_ends_is_an_error_once_the_buffer_passes_the_limit() {
    // 3 MB of one line with no newline: the buffer holds it all, because a record is only complete at its fourth newline.
    let mut chunker = FastqChunker::with_limit(1 << 20);
    let piece = vec![b'A'; 1 << 18];
    chunker.push(b"@r\n").unwrap();
    let mut failed = None;
    for i in 0..12 {
        if let Err(e) = chunker.push(&piece) {
            failed = Some((i, e));
            break;
        }
    }
    let (i, error) = failed.expect("a record of 3 MB must exceed a limit of 1 MiB");
    assert!(i <= 4, "the error comes soon after the limit, not at the end of the input (chunk {i})");
    assert!(error.to_string().contains("larger than") && error.to_string().contains("FASTQ"), "{error}");
}

#[test]
fn many_small_records_never_reach_the_limit_and_a_long_read_below_it_is_fine() {
    let mut chunker = FastqChunker::with_limit(1 << 20);
    let long = "A".repeat(500_000); // a Nanopore-sized read
    let record = format!("@long\n{long}\n+\n{}\n", "I".repeat(500_000));
    let mut n = chunker.push(record.as_bytes()).unwrap().len();
    for i in 0..20_000 {
        n += chunker.push(format!("@r{i}\nACGT\n+\nIIII\n").as_bytes()).unwrap().len();
    }
    n += chunker.finish().unwrap().len();
    assert_eq!(n, 20_001);
}

#[test]
fn the_default_limit_admits_a_record_of_a_few_megabases() {
    let mut chunker = FastqChunker::new();
    let big = format!("@big\n{}\n+\n{}\n", "A".repeat(3_000_000), "I".repeat(3_000_000));
    assert_eq!(chunker.push(big.as_bytes()).unwrap().len(), 1);
}

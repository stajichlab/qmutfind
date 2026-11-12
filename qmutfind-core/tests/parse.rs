//! FASTA and FASTQ parsing over a byte reader. No files are needed.

use qmutfind_core::{read_fasta, FastqReader};

#[test]
fn fasta_names_are_the_first_header_word_and_bases_are_uppercase_and_joined() {
    let text = ">chrA some description\nACGTac\ngtNN\n>chrB\nTTTT\n";

    let contigs = read_fasta(text.as_bytes()).unwrap();

    assert_eq!(
        contigs,
        vec![
            ("chrA".to_string(), b"ACGTACGTNN".to_vec()),
            ("chrB".to_string(), b"TTTT".to_vec()),
        ]
    );
}

#[test]
fn fastq_reader_gives_each_read_name_and_sequence_in_order() {
    let text = "@r1 extra\nACGT\n+\nIIII\n@r2\nGGGGCC\n+\nIIIIII\n";

    let reads: Vec<_> = FastqReader::new(text.as_bytes()).map(|r| r.unwrap()).collect();

    assert_eq!(
        reads,
        vec![
            ("r1".to_string(), b"ACGT".to_vec()),
            ("r2".to_string(), b"GGGGCC".to_vec()),
        ]
    );
}

#[test]
fn a_truncated_fastq_record_is_an_error_not_a_silent_drop() {
    let text = "@r1\nACGT\n+\nIIII\n@r2\nGGGG\n";

    let results: Vec<_> = FastqReader::new(text.as_bytes()).collect();

    assert!(results[0].is_ok());
    assert!(results[1].is_err());
}

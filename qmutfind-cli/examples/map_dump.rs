//! Diagnostic: map each read of a FASTQ file on its own with the qmutfind aligner and print one line per read:
//! name, mate label, number of primary alignments, then for the first primary alignment contig, 1-based start, strand, CIGAR, mapq.
//! An unmapped read has `*`. Usage: map_dump REF.fa[.gz] READS.fastq[.gz] MATE_LABEL > out.tsv
use std::io::{BufRead, BufReader, Read, Write};

use qmutfind_core::{read_fasta, Aligner, FastqChunker, ReadKind, RammapAligner};

fn open(path: &str) -> Box<dyn Read> {
    let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    if path.ends_with(".gz") { Box::new(flate2::read::MultiGzDecoder::new(f)) } else { Box::new(f) }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut fasta = Vec::new();
    open(&a[1]).read_to_end(&mut fasta).unwrap();
    let reference = read_fasta(&fasta[..]).unwrap();
    let aligner = RammapAligner::new(reference, ReadKind::Short);
    let mut reader = BufReader::new(open(&a[2]));
    let mut chunker = FastqChunker::new();
    let out = std::io::stdout();
    let mut out = std::io::BufWriter::new(out.lock());
    let mut emit = |records: Vec<(String, Vec<u8>)>| {
        for (name, seq) in records {
            let al = aligner.align(&name, &seq);
            match al.first() {
                Some(x) => writeln!(out, "{name}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", a[3], al.len(), x.contig, x.ref_start + 1, if x.reverse { '-' } else { '+' }, x.cigar_string(), x.mapq()).unwrap(),
                None => writeln!(out, "{name}\t{}\t0\t*\t0\t*\t*\t0", a[3]).unwrap(),
            }
        }
    };
    loop {
        let buf = reader.fill_buf().unwrap().to_vec();
        if buf.is_empty() { break; }
        reader.consume(buf.len());
        emit(chunker.push(&buf).unwrap());
    }
    emit(chunker.finish().unwrap());
}

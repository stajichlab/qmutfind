//! FASTA and FASTQ parsing over a byte reader. No file access and no threads,
//! so the same code serves the CLI and the browser build.

use std::io::{self, BufRead};

/// Read all FASTA records. The name is the first word of the header. Bases are uppercased.
pub fn read_fasta<R: BufRead>(reader: R) -> io::Result<Vec<(String, Vec<u8>)>> {
    let mut contigs: Vec<(String, Vec<u8>)> = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            let name = header.split_whitespace().next().unwrap_or("").to_string();
            contigs.push((name, Vec::new()));
        } else if let Some((_, seq)) = contigs.last_mut() {
            seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    Ok(contigs)
}

/// Iterator over FASTQ records as (name, sequence). The name is the first word of the header.
pub struct FastqReader<R: BufRead> {
    lines: io::Lines<R>,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        FastqReader { lines: reader.lines() }
    }

    fn next_line(&mut self, what: &str) -> io::Result<String> {
        match self.lines.next() {
            Some(line) => line,
            None => Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("truncated FASTQ record: missing {what}"),
            )),
        }
    }
}

impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = io::Result<(String, Vec<u8>)>;

    fn next(&mut self) -> Option<Self::Item> {
        let header = match self.lines.next()? {
            Ok(h) => h,
            Err(e) => return Some(Err(e)),
        };
        let record = (|| {
            let name = header
                .strip_prefix('@')
                .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "FASTQ header does not start with @"))?
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            let seq = self.next_line("sequence")?;
            self.next_line("+ line")?;
            self.next_line("quality")?;
            Ok((name, seq.trim().bytes().map(|b| b.to_ascii_uppercase()).collect()))
        })();
        Some(record)
    }
}

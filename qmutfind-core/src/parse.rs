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

/// Parse FASTQ from byte chunks that can split a record anywhere. Call `push` for each
/// chunk, then `finish` once. A truncated last record is an error.
#[derive(Default)]
pub struct FastqChunker {
    buffer: Vec<u8>,
}

impl FastqChunker {
    pub fn new() -> Self {
        FastqChunker::default()
    }

    /// Add a chunk. Returns the records that are complete now.
    pub fn push(&mut self, chunk: &[u8]) -> io::Result<Vec<(String, Vec<u8>)>> {
        self.buffer.extend_from_slice(chunk);
        // A record ends at every fourth newline.
        let mut newlines = 0;
        let mut consumed = 0;
        for (i, &b) in self.buffer.iter().enumerate() {
            if b == b'\n' {
                newlines += 1;
                if newlines % 4 == 0 {
                    consumed = i + 1;
                }
            }
        }
        let records = FastqReader::new(&self.buffer[..consumed]).collect::<io::Result<Vec<_>>>()?;
        self.buffer.drain(..consumed);
        Ok(records)
    }

    /// Finish the input. Returns a last record that had no trailing newline, or an error when
    /// the input stops inside a record.
    pub fn finish(&mut self) -> io::Result<Vec<(String, Vec<u8>)>> {
        if self.buffer.iter().all(|b| b.is_ascii_whitespace()) {
            self.buffer.clear();
            return Ok(Vec::new());
        }
        self.buffer.push(b'\n');
        let records = FastqReader::new(&self.buffer[..]).collect::<io::Result<Vec<_>>>()?;
        self.buffer.clear();
        Ok(records)
    }
}

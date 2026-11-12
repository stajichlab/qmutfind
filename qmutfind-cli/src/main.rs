//! qmutfind command-line tool.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read};
use std::process::ExitCode;

use qmutfind_core::{
    read_fasta, write_vcf, Aligner, Caller, CallerConfig, FastqReader, ReadKind, RammapAligner,
};

// Provisional defaults. They are `UNMEASURED` and are set from the benchmark (PLAN.md section 5).
const DEFAULT_MIN_DEPTH: u32 = 5;
const DEFAULT_MIN_ALT_FRACTION: f64 = 0.5;
const DEFAULT_LARGE_INDEL_MIN: usize = 50;

/// Open a plain, `.gz` or `.zst` file as a buffered reader.
fn open_input(path: &str) -> Result<Box<dyn BufRead>, String> {
    let file = File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?;
    let reader: Box<dyn Read> = if path.ends_with(".gz") {
        Box::new(flate2::read::MultiGzDecoder::new(file))
    } else if path.ends_with(".zst") {
        Box::new(zstd::stream::read::Decoder::new(file).map_err(|e| format!("{path}: {e}"))?)
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(reader)))
}

fn parse_flags(args: &[String]) -> Result<HashMap<String, String>, String> {
    let mut flags = HashMap::new();
    let mut it = args.iter();
    while let Some(key) = it.next() {
        let name = key.strip_prefix("--").ok_or_else(|| format!("unexpected argument: {key}"))?;
        let value = it.next().ok_or_else(|| format!("missing value for --{name}"))?;
        flags.insert(name.to_string(), value.clone());
    }
    Ok(flags)
}

fn required<'a>(flags: &'a HashMap<String, String>, name: &str) -> Result<&'a String, String> {
    flags.get(name).ok_or_else(|| format!("missing required option --{name}"))
}

fn call(args: &[String]) -> Result<(), String> {
    let flags = parse_flags(args)?;
    let reference_path = required(&flags, "reference")?;
    let reads_path = required(&flags, "reads")?;
    let sample = required(&flags, "sample")?;
    let out_path = required(&flags, "out")?;
    let kind = match required(&flags, "preset")?.as_str() {
        "short" => ReadKind::Short,
        "long" => ReadKind::Long,
        other => return Err(format!("--preset must be short or long, not {other}")),
    };

    let reference = read_fasta(open_input(reference_path)?).map_err(|e| format!("{reference_path}: {e}"))?;
    let contigs: Vec<(String, usize)> = reference.iter().map(|(n, s)| (n.clone(), s.len())).collect();

    let aligner = RammapAligner::new(reference.clone(), kind);
    let config = CallerConfig { min_depth: DEFAULT_MIN_DEPTH, min_alt_fraction: DEFAULT_MIN_ALT_FRACTION };
    let mut caller = Caller::new(reference, config);
    for record in FastqReader::new(open_input(reads_path)?) {
        let (name, seq) = record.map_err(|e| format!("{reads_path}: {e}"))?;
        for alignment in aligner.align(&name, &seq) {
            caller.add(&alignment, &seq);
        }
    }

    let out = File::create(out_path).map_err(|e| format!("cannot create {out_path}: {e}"))?;
    write_vcf(&mut BufWriter::new(out), &contigs, sample, &caller.calls(), DEFAULT_LARGE_INDEL_MIN)
        .map_err(|e| format!("{out_path}: {e}"))
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("call") => call(&args[1..]),
        Some(other) => Err(format!("unknown command: {other}")),
        None => Err("usage: qmutfind call --reference FA --reads FQ --sample NAME --preset short|long --out VCF".to_string()),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("qmutfind: {message}");
            ExitCode::FAILURE
        }
    }
}

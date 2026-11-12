//! qmutfind command-line tool.

mod external;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::process::ExitCode;

use qmutfind_core::{
    annotate_vcf, calls_from_freebayes_vcf, compute_stats, left_align_calls, Annotation, Genome, read_fasta, write_vcf, Caller, CallerConfig, Engine, ReadKind,
};

const DEFAULT_LARGE_INDEL_MIN: usize = 50;

/// Read a fraction flag. It must be a number from 0 to 1.
fn fraction_flag(flags: &HashMap<String, String>, name: &str, default: f64) -> Result<f64, String> {
    match flags.get(name) {
        None => Ok(default),
        Some(v) => v
            .parse::<f64>()
            .ok()
            .filter(|x| (0.0..=1.0).contains(x))
            .ok_or_else(|| format!("--{name} must be a number from 0 to 1, not {v}")),
    }
}

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

    let aligner_name = flags.get("aligner").map(String::as_str).unwrap_or("rammap");
    let caller_name = flags.get("caller").map(String::as_str).unwrap_or("builtin");
    if !matches!(aligner_name, "rammap" | "bwa-mem2") {
        return Err(format!("--aligner must be rammap or bwa-mem2, not {aligner_name}"));
    }
    if !matches!(caller_name, "builtin" | "freebayes") {
        return Err(format!("--caller must be builtin or freebayes, not {caller_name}"));
    }
    if caller_name == "freebayes" && aligner_name == "rammap" {
        return Err("freebayes needs an external aligner that gives a BAM: use --aligner bwa-mem2".to_string());
    }
    let threads: u32 = match flags.get("threads") {
        Some(v) => v.parse().map_err(|_| format!("--threads must be a whole number, not {v}"))?,
        None => 1,
    };

    let reference = read_fasta(open_input(reference_path)?).map_err(|e| format!("{reference_path}: {e}"))?;
    let contigs: Vec<(String, usize)> = reference.iter().map(|(n, s)| (n.clone(), s.len())).collect();

    let defaults = CallerConfig::defaults_for(kind);
    let config = CallerConfig {
        min_depth: match flags.get("min-depth") {
            Some(v) => v.parse().map_err(|_| format!("--min-depth must be a whole number, not {v}"))?,
            None => defaults.min_depth,
        },
        min_snp_fraction: fraction_flag(&flags, "min-snp-fraction", defaults.min_snp_fraction)?,
        min_indel_fraction: fraction_flag(&flags, "min-indel-fraction", defaults.min_indel_fraction)?,
        large_indel_min: defaults.large_indel_min,
    };
    let mut read_paths = vec![reads_path.clone()];
    if let Some(second) = flags.get("reads2") {
        read_paths.push(second.clone());
    }
    let mut caller = Caller::new(reference.clone(), config);
    let mut freebayes_calls = None;
    let mut engine_vcf: Option<String> = None;
    if aligner_name == "bwa-mem2" {
        let work = external::Workdir::create(&external::default_tmpdir(flags.get("tmpdir")))?;
        let plain_reference = external::write_plain_copy(&work, "ref.fa", &mut *open_input(reference_path)?)?;
        external::bwa_mem2_index(&work, &plain_reference)?;
        // bwa-mem2 cannot read .zst, so those files are decompressed into the work directory first.
        let mut reads: Vec<String> = Vec::new();
        for (i, path) in read_paths.iter().enumerate() {
            if path.ends_with(".zst") {
                let plain = external::write_plain_copy(&work, &format!("reads{}.fq", i + 1), &mut *open_input(path)?)?;
                reads.push(plain.to_string_lossy().into_owned());
            } else {
                reads.push(path.clone());
            }
        }
        if caller_name == "freebayes" {
            let bam = external::map_to_sorted_bam(&work, &plain_reference, &reads, sample, threads)?;
            let vcf = external::run_freebayes(&work, &plain_reference, &bam)?;
            let file = File::open(&vcf).map_err(|e| format!("cannot read {}: {e}", vcf.display()))?;
            let calls = calls_from_freebayes_vcf(BufReader::new(file), &config)
                .map_err(|e| format!("freebayes output: {e}"))?;
            // FreeBayes does not always write indels in their leftmost form.
            freebayes_calls = Some(left_align_calls(calls, &reference));
        } else {
            let mut child = external::spawn_bwa_mem2(&work, &plain_reference, &reads, sample, threads)?;
            external::for_each_sam_record(&work, &mut child, |record| caller.add_sam(record))?;
        }
    } else {
        // The rammap path runs the same Engine as the browser build, so both give the same VCF.
        let mut fasta = Vec::new();
        open_input(reference_path)?
            .read_to_end(&mut fasta)
            .map_err(|e| format!("{reference_path}: {e}"))?;
        let mut engine = Engine::new(&fasta, kind, sample)?;
        engine.set_thresholds(config.min_depth, config.min_snp_fraction, config.min_indel_fraction);
        // Version 1 maps every read, including a second mate file, as an independent single-end read.
        let mut buffer = vec![0u8; 4 << 20];
        for path in &read_paths {
            let mut reader = open_input(path)?;
            loop {
                let n = reader.read(&mut buffer).map_err(|e| format!("{path}: {e}"))?;
                if n == 0 {
                    break;
                }
                engine.add_reads(&buffer[..n]).map_err(|e| format!("{path}: {e}"))?;
            }
            engine.end_file().map_err(|e| format!("{path}: {e}"))?;
        }
        engine_vcf = Some(engine.vcf()?);
    }

    let out = File::create(out_path).map_err(|e| format!("cannot create {out_path}: {e}"))?;
    if let Some(vcf) = engine_vcf {
        return BufWriter::new(out).write_all(vcf.as_bytes()).map_err(|e| format!("{out_path}: {e}"));
    }
    let calls = freebayes_calls.unwrap_or_else(|| caller.calls());
    write_vcf(&mut BufWriter::new(out), &contigs, sample, &calls, DEFAULT_LARGE_INDEL_MIN)
        .map_err(|e| format!("{out_path}: {e}"))
}

/// `qmutfind stats`: summary statistics of a VCF, written as JSON and TSV.
fn stats(args: &[String]) -> Result<(), String> {
    let flags = parse_flags(args)?;
    let vcf_path = required(&flags, "vcf")?;
    let prefix = required(&flags, "out-prefix")?;
    let window: u64 = match flags.get("window") {
        Some(v) => v.parse().map_err(|_| format!("--window must be a whole number, not {v}"))?,
        None => 100,
    };
    let stats = compute_stats(open_input(vcf_path)?, window).map_err(|e| format!("{vcf_path}: {e}"))?;
    std::fs::write(format!("{prefix}.stats.json"), stats.to_json()).map_err(|e| format!("{prefix}.stats.json: {e}"))?;
    std::fs::write(format!("{prefix}.stats.tsv"), stats.to_tsv()).map_err(|e| format!("{prefix}.stats.tsv: {e}"))
}

/// `qmutfind annotate`: add the effect of each variant on the canonical transcript to a VCF.
fn annotate(args: &[String]) -> Result<(), String> {
    let flags = parse_flags(args)?;
    let vcf_path = required(&flags, "vcf")?;
    let reference_path = required(&flags, "reference")?;
    let out_path = required(&flags, "out")?;

    let genome = Genome::new(read_fasta(open_input(reference_path)?).map_err(|e| format!("{reference_path}: {e}"))?);
    let annotation = match flags.get("gff") {
        Some(path) => Some(Annotation::from_gff(open_input(path)?).map_err(|e| format!("{path}: {e}"))?),
        None => {
            eprintln!("qmutfind: no GFF given, so every effect is unannotated");
            None
        }
    };
    let out = File::create(out_path).map_err(|e| format!("cannot create {out_path}: {e}"))?;
    annotate_vcf(open_input(vcf_path)?, &mut BufWriter::new(out), annotation.as_ref(), &genome)
        .map_err(|e| format!("{vcf_path}: {e}"))
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("call") => call(&args[1..]),
        Some("annotate") => annotate(&args[1..]),
        Some("stats") => stats(&args[1..]),
        Some(other) => Err(format!("unknown command: {other}")),
        None => Err("usage: qmutfind call --reference FA --reads FQ --sample NAME --preset short|long --out VCF\n       qmutfind annotate --vcf VCF --reference FA [--gff GFF3] --out VCF\n       qmutfind stats --vcf VCF --out-prefix PREFIX [--window BASES]".to_string()),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("qmutfind: {message}");
            ExitCode::FAILURE
        }
    }
}

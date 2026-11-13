//! qmutfind command-line tool.

mod external;
mod fetch;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::process::ExitCode;

use qmutfind_core::{
    annotate_vcf, calls_from_freebayes_vcf, combine, CombineInput, CombineOptions, depth_plot_svg, genome_plot_svg, compute_stats, report_page, left_align_calls, Annotation, Genome, read_fasta, write_vcf, Caller, CallerConfig, Engine, ReadKind,
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
        if name == "paired" || name == "html" || name == "always-exit-zero" {
            flags.insert(name.to_string(), "true".to_string()); // a switch without a value
            continue;
        }
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

    let min_mapq: u8 = match flags.get("min-mapq") {
        Some(v) => v.parse().map_err(|_| format!("--min-mapq must be a whole number from 0 to 255, not {v}"))?,
        None => CallerConfig::default_min_mapq(kind),
    };
    let min_indel_alt_reads: u32 = match flags.get("min-indel-alt-reads") {
        Some(v) => v.parse().map_err(|_| format!("--min-indel-alt-reads must be a whole number, not {v}"))?,
        None => CallerConfig::default_min_indel_alt_reads(kind),
    };
    let float_flag = |name: &str| -> Result<f64, String> {
        match flags.get(name) {
            None => Ok(0.0),
            Some(v) => v.parse::<f64>().ok().filter(|x| x.is_finite() && *x >= 0.0).ok_or_else(|| format!("--{name} must be a number of 0 or more, not {v}")),
        }
    };
    let min_qual = float_flag("min-qual")?;
    let max_depth_ratio = float_flag("max-depth-ratio")?;
    let depth_out = flags.get("depth-out");
    let depth_window: usize = match flags.get("depth-window") {
        Some(v) => v.parse().ok().filter(|w| *w > 0).ok_or_else(|| format!("--depth-window must be a whole number of bases above 0, not {v}"))?,
        None => 10_000,
    };
    let evidence_out = flags.get("evidence-out");
    if evidence_out.is_some() && caller_name == "freebayes" {
        return Err("--evidence-out needs the built-in caller: the FreeBayes path keeps no per-base counts".to_string());
    }
    if depth_out.is_some() && caller_name == "freebayes" {
        return Err("--depth-out needs the built-in caller: the FreeBayes path keeps no per-base counts".to_string());
    }
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
    let paired = flags.contains_key("paired");
    if paired && !flags.contains_key("reads2") {
        return Err("--paired needs --reads2: the second file of the pairs".to_string());
    }
    let max_insert: usize = match flags.get("max-insert") {
        Some(v) => v.parse().ok().filter(|m| *m > 0).ok_or_else(|| format!("--max-insert must be a whole number of bases above 0, not {v}"))?,
        None => qmutfind_core::DEFAULT_MAX_INSERT,
    };
    let mut caller = Caller::new(reference.clone(), config);
    caller.set_min_mapq(min_mapq);
    caller.set_min_indel_alt_reads(min_indel_alt_reads);
    caller.set_min_qual(min_qual);
    caller.set_max_depth_ratio(max_depth_ratio);
    let mut freebayes_calls = None;
    let mut engine_vcf: Option<String> = None;
    let mut depth_bed: Option<String> = None;
    let mut evidence_text: Option<String> = None;
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
            let vcf = external::run_freebayes(&work, &plain_reference, &bam, min_mapq)?;
            let file = File::open(&vcf).map_err(|e| format!("cannot read {}: {e}", vcf.display()))?;
            let calls = calls_from_freebayes_vcf(BufReader::new(file), &config)
                .map_err(|e| format!("freebayes output: {e}"))?;
            // FreeBayes does not always write indels in their leftmost form.
            let calls: Vec<_> = calls
                .into_iter()
                .filter(|c| c.reference.len() == c.alt.len() || c.alt_count >= min_indel_alt_reads)
                .collect();
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
        engine.set_min_mapq(min_mapq);
        engine.set_min_indel_alt_reads(min_indel_alt_reads);
        engine.set_min_qual(min_qual);
        engine.set_max_depth_ratio(max_depth_ratio);
        let mut buffer = vec![0u8; 4 << 20];
        if paired {
            // Both files are read in turn, one chunk each, and mapped as pairs (see `qmutfind_core::choose_pair`).
            engine.set_max_insert(max_insert);
            let mut readers = [open_input(&read_paths[0])?, open_input(&read_paths[1])?];
            let mut done = [false, false];
            while !(done[0] && done[1]) {
                for k in 0..2 {
                    if done[k] {
                        continue;
                    }
                    let n = readers[k].read(&mut buffer).map_err(|e| format!("{}: {e}", read_paths[k]))?;
                    if n == 0 {
                        done[k] = true;
                    } else if k == 0 {
                        engine.add_mate1(&buffer[..n]).map_err(|e| format!("{}: {e}", read_paths[k]))?;
                    } else {
                        engine.add_mate2(&buffer[..n]).map_err(|e| format!("{}: {e}", read_paths[k]))?;
                    }
                }
            }
            engine.end_mates().map_err(|e| format!("{} and {}: {e}", read_paths[0], read_paths[1]))?;
        }
        // Without --paired every read, including a second file, is mapped as an independent single-end read.
        for path in read_paths.iter().filter(|_| !paired) {
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
        depth_bed = depth_out.map(|_| engine.depth_bed(depth_window));
        evidence_text = evidence_out.map(|_| engine.evidence_vcf()).transpose()?;
        engine_vcf = Some(engine.vcf()?);
    }

    if let Some(path) = depth_out {
        let bed = depth_bed.unwrap_or_else(|| caller.depth_bed(depth_window));
        std::fs::write(path, bed).map_err(|e| format!("{path}: {e}"))?;
    }
    if let Some(path) = evidence_out {
        let text = match evidence_text {
            Some(text) => text,
            None => {
                let mut buffer = Vec::new();
                write_vcf(&mut buffer, &contigs, sample, &caller.evidence_calls(), DEFAULT_LARGE_INDEL_MIN).map_err(|e| format!("{path}: {e}"))?;
                String::from_utf8(buffer).map_err(|e| e.to_string())?
            }
        };
        std::fs::write(path, text).map_err(|e| format!("{path}: {e}"))?;
    }
    let out = File::create(out_path).map_err(|e| format!("cannot create {out_path}: {e}"))?;
    if let Some(vcf) = engine_vcf {
        return BufWriter::new(out).write_all(vcf.as_bytes()).map_err(|e| format!("{out_path}: {e}"));
    }
    let calls = freebayes_calls.unwrap_or_else(|| caller.calls());
    write_vcf(&mut BufWriter::new(out), &contigs, sample, &calls, DEFAULT_LARGE_INDEL_MIN)
        .map_err(|e| format!("{out_path}: {e}"))
}

/// `qmutfind report`: a page with the gene models and the variants of one region.
fn report(args: &[String]) -> Result<(), String> {
    let flags = parse_flags(args)?;
    let vcf_path = required(&flags, "vcf")?;
    let region = required(&flags, "region")?;
    let out_path = required(&flags, "out")?;
    let title = flags.get("title").cloned().unwrap_or_else(|| format!("Mutations in {region}"));
    let mut vcf = String::new();
    open_input(vcf_path)?.read_to_string(&mut vcf).map_err(|e| format!("{vcf_path}: {e}"))?;
    let gff = match flags.get("gff") {
        Some(path) => {
            let mut bytes = Vec::new();
            open_input(path)?.read_to_end(&mut bytes).map_err(|e| format!("{path}: {e}"))?;
            Some(bytes)
        }
        None => None,
    };
    let page = report_page(&vcf, gff.as_deref(), region, &title)?;
    std::fs::write(out_path, page).map_err(|e| format!("{out_path}: {e}"))
}

/// `qmutfind plot`: the whole-genome plot of a VCF as an SVG file. Contig lengths come from the VCF's
/// `##contig` lines, or from `--reference` when it is given.
fn plot(args: &[String]) -> Result<(), String> {
    let flags = parse_flags(args)?;
    let vcf_path = required(&flags, "vcf")?;
    let out_path = required(&flags, "out")?;
    let mut vcf = String::new();
    open_input(vcf_path)?.read_to_string(&mut vcf).map_err(|e| format!("{vcf_path}: {e}"))?;
    let contigs = match flags.get("reference") {
        Some(path) => {
            let reference = read_fasta(open_input(path)?).map_err(|e| format!("{path}: {e}"))?;
            Some(reference.into_iter().map(|(name, seq)| (name, seq.len())).collect::<Vec<_>>())
        }
        None => None,
    };
    let svg = genome_plot_svg(&vcf, contigs.as_deref(), flags.get("order").map(String::as_str).unwrap_or("input"))?;
    std::fs::write(out_path, svg).map_err(|e| format!("{out_path}: {e}"))
}

/// `qmutfind depth-plot`: the genome-wide depth plot from BED text (our `--depth-out` file or `mosdepth` output).
fn depth_plot(args: &[String]) -> Result<(), String> {
    let flags = parse_flags(args)?;
    let bed_path = required(&flags, "bed")?;
    let out_path = required(&flags, "out")?;
    let mut bed = String::new();
    open_input(bed_path)?.read_to_string(&mut bed).map_err(|e| format!("{bed_path}: {e}"))?;
    let vcf = match flags.get("vcf") {
        Some(path) => {
            let mut text = String::new();
            open_input(path)?.read_to_string(&mut text).map_err(|e| format!("{path}: {e}"))?;
            Some(text)
        }
        None => None,
    };
    let contigs = match flags.get("reference") {
        Some(path) => {
            let reference = read_fasta(open_input(path)?).map_err(|e| format!("{path}: {e}"))?;
            Some(reference.into_iter().map(|(name, seq)| (name, seq.len())).collect::<Vec<_>>())
        }
        None => None,
    };
    let order = flags.get("order").map(String::as_str).unwrap_or("input");
    let svg = depth_plot_svg(&bed, contigs.as_deref(), order, vcf.as_deref())?;
    std::fs::write(out_path, svg).map_err(|e| format!("{out_path}: {e}"))
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

/// `qmutfind combine`: one table for several strains, with an optional parent. `--strain NAME=VCF[,BED]` is repeated.
fn combine_command(args: &[String]) -> Result<(), String> {
    let mut strains: Vec<(String, String, Option<String>, Option<String>)> = Vec::new();
    let mut parent: Option<String> = None;
    let mut prefix: Option<String> = None;
    let mut min_depth = 5.0; // the caller's own minimum depth for short reads; the coverage windows are coarse
    let mut it = args.iter();
    while let Some(key) = it.next() {
        let name = key.strip_prefix("--").ok_or_else(|| format!("unexpected argument: {key}"))?;
        let value = it.next().ok_or_else(|| format!("missing value for --{name}"))?;
        match name {
            "strain" => {
                let (strain, files) = value.split_once('=').filter(|(n, f)| !n.is_empty() && !f.is_empty()).ok_or_else(|| format!("--strain must look like NAME=VCF, NAME=VCF,BED or NAME=VCF,BED,EVIDENCE_VCF, not {value}"))?;
                let parts: Vec<&str> = files.split(',').collect();
                if parts.len() > 3 || parts[0].is_empty() {
                    return Err(format!("--strain must look like NAME=VCF, NAME=VCF,BED or NAME=VCF,BED,EVIDENCE_VCF, not {value}"));
                }
                let optional = |i: usize| parts.get(i).filter(|p| !p.is_empty()).map(|p| p.to_string());
                strains.push((strain.to_string(), parts[0].to_string(), optional(1), optional(2)));
            }
            "parent" => parent = Some(value.clone()),
            "out-prefix" => prefix = Some(value.clone()),
            "min-depth" => min_depth = value.parse().map_err(|_| format!("--min-depth must be a number, not {value}"))?,
            other => return Err(format!("unknown option --{other}")),
        }
    }
    if strains.is_empty() {
        return Err("give at least one --strain NAME=VCF[,BED]".to_string());
    }
    let prefix = prefix.ok_or("missing required option --out-prefix")?;
    let read = |path: &str| -> Result<String, String> {
        let mut text = String::new();
        open_input(path)?.read_to_string(&mut text).map_err(|e| format!("{path}: {e}"))?;
        Ok(text)
    };
    let mut texts = Vec::new();
    for (name, vcf, bed, evidence) in &strains {
        texts.push((name.as_str(), read(vcf)?, bed.as_deref().map(read).transpose()?, evidence.as_deref().map(read).transpose()?));
    }
    let inputs: Vec<CombineInput> = texts
        .iter()
        .map(|(name, vcf, bed, evidence)| CombineInput { name, vcf, depth_bed: bed.as_deref(), evidence_vcf: evidence.as_deref() })
        .collect();
    let combined = combine(&inputs, &CombineOptions { parent: parent.as_deref(), min_depth })?;
    for warning in &combined.warnings {
        eprintln!("qmutfind: warning: {warning}");
    }
    std::fs::write(format!("{prefix}.combined.tsv"), combined.to_tsv()).map_err(|e| format!("{prefix}.combined.tsv: {e}"))?;
    std::fs::write(format!("{prefix}.summary.tsv"), combined.summary_tsv()).map_err(|e| format!("{prefix}.summary.tsv: {e}"))?;
    std::fs::write(format!("{prefix}.combined.vcf"), combined.to_vcf()).map_err(|e| format!("{prefix}.combined.vcf: {e}"))
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

/// `qmutfind context --reference FA [--gff GFF3] --variant CONTIG:POS:REF:ALT [--flank N] [--html] [--out FILE]`: the sequence around a variant.
fn context(args: &[String]) -> Result<(), String> {
    let flags = parse_flags(args)?;
    let reference_path = required(&flags, "reference")?;
    let variant = required(&flags, "variant")?;
    let parts: Vec<&str> = variant.split(':').collect();
    let [contig, pos, reference, alt] = parts[..] else {
        return Err(format!("--variant must be CONTIG:POS:REF:ALT (for example chr1:50:A:G), not {variant}"));
    };
    let pos: u64 = pos.parse().map_err(|_| format!("the position in --variant must be a whole number, not {pos}"))?;
    let flank: u64 = match flags.get("flank") {
        Some(v) => v.parse().map_err(|_| format!("--flank must be a whole number of bases, not {v}"))?,
        None => 40,
    };
    let mut fasta = Vec::new();
    open_input(reference_path)?.read_to_end(&mut fasta).map_err(|e| format!("{reference_path}: {e}"))?;
    let gff = match flags.get("gff") {
        Some(path) => {
            let mut bytes = Vec::new();
            open_input(path)?.read_to_end(&mut bytes).map_err(|e| format!("{path}: {e}"))?;
            Some(bytes)
        }
        None => None,
    };
    let text = qmutfind_core::context_text(&fasta, gff.as_deref(), contig, pos, reference, alt, flank, flags.contains_key("html"))?;
    match flags.get("out") {
        Some(path) => std::fs::write(path, text).map_err(|e| format!("{path}: {e}")),
        None => {
            print!("{text}");
            Ok(())
        }
    }
}

/// `qmutfind confirm`: is an introduced cassette integrated where it should be, and intact? Reads are mapped to the genome plus the cassette;
/// see `qmutfind_core::Confirmer`. Exit status: 0 confirmed, 3 not confirmed, 4 inconclusive (1 for an error); `--always-exit-zero` gives 0 after a report.
fn confirm_command(args: &[String]) -> Result<(), String> {
    use qmutfind_core::{find_gene, ConfirmOptions, Confirmer, Target, Verdict};
    let flags = parse_flags(args)?;
    let reference_path = required(&flags, "reference")?;
    let construct_path = flags.get("construct").ok_or_else(|| "missing required option --construct: a FASTA file with the sequence you introduced (for example the hph cassette)".to_string())?;
    let reads_path = required(&flags, "reads")?;
    let prefix = required(&flags, "out-prefix")?;
    let kind = match required(&flags, "preset")?.as_str() {
        "short" => ReadKind::Short,
        "long" => ReadKind::Long,
        other => return Err(format!("--preset must be short or long, not {other}")),
    };
    let number = |name: &str, default: f64| -> Result<f64, String> {
        match flags.get(name) {
            None => Ok(default),
            Some(v) => v.parse::<f64>().ok().filter(|x| *x >= 0.0).ok_or_else(|| format!("--{name} must be a number of 0 or more, not {v}")),
        }
    };
    let mut genome = read_fasta(open_input(reference_path)?).map_err(|e| format!("{reference_path}: {e}"))?;
    let construct = read_fasta(open_input(construct_path)?).map_err(|e| format!("{construct_path}: {e}"))?;
    if construct.is_empty() {
        return Err(format!("{construct_path} has no sequence"));
    }
    for (name, _) in &construct {
        if genome.iter().any(|(g, _)| g == name) {
            return Err(format!("the construct contig {name} has the same name as a contig of the reference; rename it in {construct_path}"));
        }
    }
    let names: Vec<String> = construct.iter().map(|(n, _)| n.clone()).collect();
    genome.extend(construct);
    let mut fasta = Vec::new();
    for (name, seq) in &genome {
        fasta.extend_from_slice(format!(">{name}\n").as_bytes());
        fasta.extend_from_slice(seq);
        fasta.push(b'\n');
    }
    drop(genome);
    let target = match (flags.get("target"), flags.get("gene")) {
        (Some(_), Some(_)) => return Err("give --target or --gene, not both".to_string()),
        (Some(t), None) => {
            let bad = || format!("--target must be CONTIG:START-END (1-based, for example chr1:50001-51500), not {t}");
            let (contig, range) = t.rsplit_once(':').ok_or_else(bad)?;
            let (start, end) = range.replace(',', "").split_once('-').map(|(a, b)| (a.to_string(), b.to_string())).ok_or_else(bad)?;
            let (start, end): (u64, u64) = (start.parse().map_err(|_| bad())?, end.parse().map_err(|_| bad())?);
            if start == 0 || end < start {
                return Err(bad());
            }
            Some(Target { contig: contig.to_string(), start, end })
        }
        (None, Some(gene)) => {
            let gff = flags.get("gff").ok_or_else(|| "--gene needs --gff: the GFF3 file that holds the gene".to_string())?;
            Some(find_gene(open_input(gff)?, gene).map_err(|e| format!("{gff}: {e}"))?)
        }
        (None, None) => None,
    };
    let defaults = ConfirmOptions::default();
    let options = ConfirmOptions {
        target,
        min_junction_reads: number("min-junction-reads", defaults.min_junction_reads as f64)? as u32,
        min_depth: number("min-depth", defaults.min_depth)?,
        target_tolerance: number("tolerance", defaults.target_tolerance as f64)? as u64,
        ..defaults
    };
    let constructs: Vec<&str> = names.iter().map(String::as_str).collect();
    let mut confirmer = Confirmer::new(&fasta, &constructs, kind, options)?;
    drop(fasta);
    let mut buffer = vec![0u8; 4 << 20];
    let mut paths = vec![reads_path.clone()];
    if let Some(second) = flags.get("reads2") {
        paths.push(second.clone());
    }
    for path in &paths {
        let mut reader = open_input(path)?;
        loop {
            let n = reader.read(&mut buffer).map_err(|e| format!("{path}: {e}"))?;
            if n == 0 {
                break;
            }
            confirmer.add_reads(&buffer[..n]).map_err(|e| format!("{path}: {e}"))?;
        }
        confirmer.end_file().map_err(|e| format!("{path}: {e}"))?;
    }
    let report = confirmer.report();
    let write = |suffix: &str, text: String| std::fs::write(format!("{prefix}.{suffix}"), text).map_err(|e| format!("{prefix}.{suffix}: {e}"));
    write("confirm.txt", report.to_text())?;
    write("confirm.json", report.to_json())?;
    let side = |s: qmutfind_core::Side| if s == qmutfind_core::Side::Left { "left" } else { "right" };
    let mut tsv = String::from("contig_a\tpos_a\tside_a\tcontig_b\tpos_b\tside_b\treads\toverlap\tkind\n");
    for j in &report.junctions {
        let (ca, cb) = (names.contains(&j.a.contig), names.contains(&j.b.contig));
        let kind = match (ca, cb) { (true, true) => "cassette_internal", (false, false) => "genome", _ => "integration" };
        tsv.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{kind}\n", j.a.contig, j.a.pos, side(j.a.keeps), j.b.contig, j.b.pos, side(j.b.keeps), j.reads, j.overlap));
    }
    write("junctions.tsv", tsv)?;
    let mut depth = String::from("contig\tpos\tdepth\n");
    for k in &report.cassettes {
        for (i, d) in k.depth.iter().enumerate() {
            depth.push_str(&format!("{}\t{}\t{d}\n", k.name, i + 1));
        }
    }
    write("cassette_depth.tsv", depth)?;
    print!("{}", report.to_text());
    std::io::stdout().flush().map_err(|e| e.to_string())?;
    let code = match report.verdict() { Verdict::Confirmed => 0, Verdict::NotConfirmed => 3, Verdict::Inconclusive => 4 };
    if code != 0 && !flags.contains_key("always-exit-zero") {
        std::process::exit(code);
    }
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let result = match args.first().map(String::as_str) {
        Some("call") => call(&args[1..]),
        Some("annotate") => annotate(&args[1..]),
        Some("context") => context(&args[1..]),
        Some("confirm") => confirm_command(&args[1..]),
        Some("stats") => stats(&args[1..]),
        Some("combine") => combine_command(&args[1..]),
        Some("report") => report(&args[1..]),
        Some("plot") => plot(&args[1..]),
        Some("depth-plot") => depth_plot(&args[1..]),
        Some("fetch") => fetch::fetch(&args[1..]),
        Some(other) => Err(format!("unknown command: {other}")),
        None => Err("usage: qmutfind call --reference FA --reads FQ --sample NAME --preset short|long --out VCF\n       qmutfind annotate --vcf VCF --reference FA [--gff GFF3] --out VCF\n       qmutfind context --reference FA [--gff GFF3] --variant CONTIG:POS:REF:ALT [--flank BASES] [--html] [--out FILE]\n       qmutfind confirm --reference FA --construct FA --reads FQ [--reads2 FQ] --preset short|long (--target CONTIG:START-END | --gene NAME --gff GFF3) --out-prefix PREFIX\n       qmutfind stats --vcf VCF --out-prefix PREFIX [--window BASES]\n       qmutfind combine --strain NAME=VCF[,BED[,EVIDENCE_VCF]] ... [--parent NAME] [--min-depth DEPTH] --out-prefix PREFIX\n       qmutfind report --vcf VCF [--gff GFF3] --region CONTIG:START-END [--title TEXT] --out HTML\n       qmutfind plot --vcf VCF [--reference FA] [--order input|length|name] --out SVG\n       qmutfind depth-plot --bed BED [--vcf VCF] [--reference FA] [--order input|length|name] --out SVG\n       qmutfind fetch --organism NAME --out-dir DIR [--eupathdb-key-file FILE]   (the FASTA needs a key: EUPATHDB_KEY or --eupathdb-key-file)".to_string()),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(message) => {
            eprintln!("qmutfind: {message}");
            ExitCode::FAILURE
        }
    }
}

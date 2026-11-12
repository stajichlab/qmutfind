//! Run external programs (bwa-mem2, samtools, freebayes) for the selectable backends.
//! The programs come from PATH. The pinned versions are in the pixi environment.

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};

use qmutfind_core::{parse_sam_line, SamRecord};

/// A temporary directory that is removed when it goes out of scope, also after an error.
pub struct Workdir {
    path: PathBuf,
}

impl Workdir {
    pub fn create(parent: &Path) -> Result<Workdir, String> {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.subsec_nanos())
            .unwrap_or(0);
        let path = parent.join(format!("qmutfind-{}-{}", std::process::id(), nanos));
        std::fs::create_dir_all(&path).map_err(|e| format!("cannot create {}: {e}", path.display()))?;
        Ok(Workdir { path })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for Workdir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.path);
    }
}

/// Where temporary files go: --tmpdir, then $SCRATCH, then $TMPDIR, then the system default.
pub fn default_tmpdir(flag: Option<&String>) -> PathBuf {
    if let Some(dir) = flag {
        return PathBuf::from(dir);
    }
    for var in ["SCRATCH", "TMPDIR"] {
        if let Ok(dir) = std::env::var(var) {
            if !dir.is_empty() {
                return PathBuf::from(dir);
            }
        }
    }
    std::env::temp_dir()
}

fn spawn_error(tool: &str, e: std::io::Error) -> String {
    if e.kind() == std::io::ErrorKind::NotFound {
        format!("{tool} was not found on PATH (the pinned tools are in the pixi environment: run through `pixi run`)")
    } else {
        format!("cannot start {tool}: {e}")
    }
}

fn log_tail(log: &Path) -> String {
    let text = std::fs::read_to_string(log).unwrap_or_default();
    let lines: Vec<&str> = text.lines().collect();
    lines[lines.len().saturating_sub(8)..].join("\n")
}

/// Run a program to completion. Its standard error goes to a log file in `work`.
pub fn run_tool(work: &Workdir, tool: &str, args: &[&str]) -> Result<(), String> {
    let log = work.path().join(format!("{tool}.log"));
    let status = Command::new(tool)
        .args(args)
        .stdout(Stdio::null())
        .stderr(Stdio::from(File::create(&log).map_err(|e| e.to_string())?))
        .status()
        .map_err(|e| spawn_error(tool, e))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{tool} failed ({status}):\n{}", log_tail(&log)))
    }
}

/// Write the reference as plain FASTA into `work`. The tools cannot read `.zst`, and
/// a plain copy also serves a `.gz` input.
pub fn write_plain_copy(work: &Workdir, name: &str, reader: &mut dyn BufRead) -> Result<PathBuf, String> {
    let dest = work.path().join(name);
    let mut out = File::create(&dest).map_err(|e| format!("cannot write {}: {e}", dest.display()))?;
    std::io::copy(reader, &mut out).map_err(|e| format!("cannot copy to {}: {e}", dest.display()))?;
    out.flush().map_err(|e| e.to_string())?;
    Ok(dest)
}

/// Index the reference with bwa-mem2.
pub fn bwa_mem2_index(work: &Workdir, reference: &Path) -> Result<(), String> {
    run_tool(work, "bwa-mem2", &["index", &reference.to_string_lossy()])
}

/// Start `bwa-mem2 mem`. SAM text comes out on the child's standard output.
pub fn spawn_bwa_mem2(
    work: &Workdir,
    reference: &Path,
    reads: &[String],
    sample: &str,
    threads: u32,
) -> Result<Child, String> {
    let log = work.path().join("bwa-mem2-mem.log");
    let read_group = format!("@RG\\tID:{sample}\\tSM:{sample}");
    Command::new("bwa-mem2")
        .args(["mem", "-t", &threads.to_string(), "-R", &read_group])
        .arg(reference)
        .args(reads)
        .stdout(Stdio::piped())
        .stderr(Stdio::from(File::create(&log).map_err(|e| e.to_string())?))
        .spawn()
        .map_err(|e| spawn_error("bwa-mem2", e))
}

/// Give every mapped SAM record on the child's standard output to `f`, then check the exit status.
pub fn for_each_sam_record(
    work: &Workdir,
    child: &mut Child,
    mut f: impl FnMut(&SamRecord),
) -> Result<(), String> {
    let stdout = child.stdout.take().ok_or("no standard output from the aligner")?;
    for line in BufReader::new(stdout).lines() {
        let line = line.map_err(|e| format!("reading SAM from the aligner: {e}"))?;
        if let Some(record) = parse_sam_line(&line)? {
            f(&record);
        }
    }
    let status = child.wait().map_err(|e| e.to_string())?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("bwa-mem2 failed ({status}):\n{}", log_tail(&work.path().join("bwa-mem2-mem.log"))))
    }
}

/// Run a program to completion with its standard output written to `out`.
pub fn run_tool_to_file(work: &Workdir, tool: &str, args: &[&str], out: &Path) -> Result<(), String> {
    let log = work.path().join(format!("{tool}-out.log"));
    let status = Command::new(tool)
        .args(args)
        .stdout(Stdio::from(File::create(out).map_err(|e| format!("cannot write {}: {e}", out.display()))?))
        .stderr(Stdio::from(File::create(&log).map_err(|e| e.to_string())?))
        .status()
        .map_err(|e| spawn_error(tool, e))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{tool} failed ({status}):\n{}", log_tail(&log)))
    }
}

/// Map with bwa-mem2, sort with samtools, and index. Returns the sorted BAM.
pub fn map_to_sorted_bam(
    work: &Workdir,
    reference: &Path,
    reads: &[String],
    sample: &str,
    threads: u32,
) -> Result<PathBuf, String> {
    let mut bwa = spawn_bwa_mem2(work, reference, reads, sample, threads)?;
    let bam = work.path().join("aln.bam");
    let sort_log = work.path().join("samtools-sort.log");
    let mut sort = Command::new("samtools")
        .args(["sort", "-@", &threads.to_string(), "-o"])
        .arg(&bam)
        .arg("-")
        .stdin(Stdio::from(bwa.stdout.take().ok_or("no output from the aligner")?))
        .stdout(Stdio::null())
        .stderr(Stdio::from(File::create(&sort_log).map_err(|e| e.to_string())?))
        .spawn()
        .map_err(|e| spawn_error("samtools", e))?;
    let sort_status = sort.wait().map_err(|e| e.to_string())?;
    let bwa_status = bwa.wait().map_err(|e| e.to_string())?;
    if !bwa_status.success() {
        return Err(format!("bwa-mem2 failed ({bwa_status}):\n{}", log_tail(&work.path().join("bwa-mem2-mem.log"))));
    }
    if !sort_status.success() {
        return Err(format!("samtools sort failed ({sort_status}):\n{}", log_tail(&sort_log)));
    }
    run_tool(work, "samtools", &["index", &bam.to_string_lossy()])?;
    Ok(bam)
}

/// Run FreeBayes for a haploid sample. Returns the raw VCF.
pub fn run_freebayes(work: &Workdir, reference: &Path, bam: &Path, min_mapq: u8) -> Result<PathBuf, String> {
    run_tool(work, "samtools", &["faidx", &reference.to_string_lossy()])?;
    let vcf = work.path().join("freebayes.vcf");
    run_tool_to_file(
        work,
        "freebayes",
        &[
            "-f",
            &reference.to_string_lossy(),
            "--ploidy",
            "1",
            "--min-mapping-quality",
            &min_mapq.max(1).to_string(),
            &bam.to_string_lossy(),
        ],
        &vcf,
    )?;
    Ok(vcf)
}

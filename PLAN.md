# qmutfind: plan

Licence (user decision, 2026-11-11): GPL-3.0. `LICENSE` holds the standard text (35,149 bytes,
byte-identical to the copy in pinbo/bwa-samtools-web). The Cargo manifests declare `GPL-3.0-only`.
ASSUMED: "only" and not "or later". Change `license` in `Cargo.toml` if you want "or later".
`rammap` declares `license = "MIT"` in its workspace `Cargo.toml` (read 2026-11-11). I did not read its
`LICENSE` file. The legal effect of combining MIT and GPL-3.0 code is not checked here.

Status (2026-11-12). Works and is tested (95 Rust tests, 25 Python tests, 2 Emscripten tests, 7 Node
WASM tests; CI passed on GitHub on 2026-11-12):
- `qmutfind call` with selectable backends: `--aligner rammap|bwa-mem2`, `--caller builtin|freebayes`,
  `--reads2`, `--threads`, `--min-depth`, `--min-snp-fraction`, `--min-indel-fraction`. The `rammap` path runs the
  same `Engine` as the browser build.
- `qmutfind annotate`: the effect of each variant on the canonical transcript of a gene, for any VCF
  (`EFFECT`, `GENE`, `TRANSCRIPT`, `AA`, `NTX`). The tag names and the canonical-transcript rule are ASSUMED:
  the user did not answer question 43 before going offline.
- WebAssembly: `Session` (reference bytes, read chunks, VCF out) built with `pixi run build-wasm-node` or
  `build-wasm-web`; `pixi run test-wasm` checks that the WASM VCF equals the CLI VCF. Emscripten 3.1.64 is
  installed by `pixi run setup-emsdk`.
Not built: statistics, the report (SVG locus views), the "possible large indel, unresolved" report for short
reads, mapping-quality and base-quality filters, WASM bindings for `annotate`, the genomics-wasm BWA-MEM2 and
FreeBayes builds, a real-browser test. Defaults for thresholds were set from data and are tuned on the data
they were measured on (see section 11).
The original voice-note text is kept in `PLAN.original-voice-notes.md`.

Markers used in this file:
- `UNMEASURED`: a number or fact with no data behind it. A measurement must settle it.
- `UNVERIFIED`: a statement taken from memory or from the user. It was not checked.
- `ASSUMED`: a default I chose. The user did not state it. Change it if wrong.

## 1. Goal

Find mutations in one fungal sample by mapping its reads to a reference genome.
The tool runs in two forms from one code base:
- a standalone command-line tool (CLI);
- a browser tool built with WebAssembly (WASM).

Input reads are Illumina short reads or long reads (Nanopore, PacBio).
Reads can be compressed FASTQ.

## 2. Scope of version 1

- Each sample is haploid.
- Each sample is processed on its own, one at a time, against a reference genome.
- The tool then merges the per-sample results into one composite result.
  Every record carries a sample label.
- Parent-versus-mutant subtraction is not built in. A user can filter on `NSAMPLES`
  (section 7) to approximate it.

## 3. Architecture

A Cargo workspace with three crates.

| Crate | Role |
|---|---|
| `qmutfind-core` | Pure Rust logic. Holds the `Aligner` trait, the pileup caller, the effect classifier and the statistics. It does no file I/O and starts no threads. It compiles for `wasm32`. It takes reads as a stream, one batch at a time. |
| `qmutfind-cli` | Native binary. Reads plain and compressed files (`.gz`, `.zst`). May use threads. Subcommands: `call` and `report`. |
| `qmutfind-wasm` | Thin `wasm-bindgen` wrapper. Runs in a Web Worker. JavaScript reads the FASTQ in chunks with the File API and passes the chunks to the core. |

Rules:
- The CLI and the browser build must give the same VCF for the same small input.
  A test checks this.
- The first WASM build may be single-threaded. WASM threads need extra browser features.
  `UNMEASURED`: whether `rammap` needs threads, and how long a real *A. fumigatus*
  run takes in a browser.
- WASM32 has a 4 GB memory limit. Processing one sample at a time keeps memory bounded.
- New Python scripts that read large text files accept plain and compressed input.
  (Global lab rule. Applies to any helper script in this project.)

## 4. Aligner

The core depends on a small `Aligner` trait.
It takes reads and a reference. It returns alignments with CIGAR strings.

Decision: test `rammap` first.
- `rammap` is a pure-Rust aligner. It mirrors the `minimap2` interface and output.
  It has presets for short reads and long reads. Source: https://github.com/jwanglab/rammap
- `UNVERIFIED`: whether `rammap` builds for `wasm32-unknown-unknown`, its maturity,
  and its dependencies (threads, C code).

Spike (run 2026-11-11 on a SLURM node: 4 CPUs, 8 GB; `rammap` v1.2.0, commit 3702b34, `minimap2` 2.30,
Rust 1.96.0, Node 22.23.1). The test script is `spike/rammap/wasm_run.cjs`. Data are outside the repository
in `/bigdata/stajichlab/jstajich/projects/qmutfind_data` (Af293 FASTA and GFF3 with MD5 checked, three runs
with read counts checked against ENA).

Result: `rammap` passes every test that was run. Details:

| Test | Result |
|---|---|
| Native build (`cargo build --release -p rammap`) | OK. 55.9 s, peak 420 MB. |
| WASM build (`-p rammap-core --no-default-features --target wasm32-unknown-unknown`, `RUSTFLAGS="-C target-feature=+simd128"`) | OK. 25.9 s, peak 481 MB. `.wasm` is 861 KB raw and 425 KB after `wasm-bindgen`. |
| Native versus `minimap2`, Illumina (SRR6277229, 1,427,487 pairs, `-x sr`, 4 threads) | All 2,855,733 SAM records identical on fields 1-9 and the `NM` tag. 2,400,691 mapped (84.07%) in both. `minimap2` 14.6 s, 865 MB. `rammap` 19.4 s, 1.26 GB. |
| Native versus `minimap2`, Nanopore (SRR28035789, 624,929 reads, `-x map-ont`, 4 threads) | 739,669 records identical including `NM` and `AS` (99.97%). `minimap2` only: 210. `rammap` only: 194. Of these, 207 and 191 are secondary or supplementary records. 3 reads differ in the primary record, and each of the 3 differs in MAPQ (for example 27 versus 28). Unmapped reads are equal (63,871). `minimap2` 2:21.8, 4.09 GB. `rammap` 2:13.1, 4.11 GB. |
| WASM (Node) versus native `rammap`, Illumina, first 300k reads, single-end, 1 thread | 300,089 of 300,089 records identical on fields 1-6. WASM 6.7 s (index 1.4 s, align 4.9 s), peak 377 MB. Native 4.9 s, peak 532 MB. |
| WASM (Node) versus native `rammap`, Nanopore, first 100k reads, 1 thread | 117,836 of 117,836 records identical on fields 1-6. WASM 86.8 s (index 1.5 s), peak 555 MB. Native 46.8 s, peak 1.06 GB. |

Notes and cautions:
- The first WASM build did not use `+simd128`. The repository sets the flag in `.cargo/config.toml`, and
  Cargo reads that file from the current directory, not from `--manifest-path`. Pass `RUSTFLAGS` explicitly.
- The `wasm-bindgen` CLI version must equal the library version in `Cargo.lock` (0.2.114). The CLI
  installed on this machine was 0.2.130 and would not match.
- WASM was 1.4x slower than native single-thread on short reads and 1.9x slower on Nanopore reads
  (these are single measurements on one machine).
- Extrapolation, `UNMEASURED`: a full Nanopore run of 625k reads would take about 9 minutes in
  single-threaded WASM if time scales linearly with read count.
- Not tested: a real browser (only Node, which uses the same V8 engine), the multi-threaded build (needs
  nightly Rust and `SharedArrayBuffer`), gzip reference input, runs larger than 100k reads in WASM,
  Nanopore comparison against `minimap2` on a second genome.
- `AlignSession` (the streaming WASM API) takes one FASTQ stream. It has no paired-end mode.
  Short-read pairs would be mapped as independent single-end reads in the browser. The native CLI
  maps pairs. This is a gap for the browser build (see open items).
- `rammap` also has a Rust library API (`Aligner`, `map_seq`) that is compatible with `minimap2-rs`.
  `UNVERIFIED`: that this API builds and runs for `wasm32` without the `wasm-bindgen` entry points.
  Our `qmutfind-core` would use the library API, not `AlignSession`.

Decision consequence: keep `rammap` as the aligner. The C `minimap2` fallback is not needed unless a
browser test fails. Keep the `Aligner` trait anyway, because it costs little.

Fallback if the spike fails: C `minimap2` built with Emscripten or WASI.
This keeps the browser build but loses the "pure Rust, one build" goal.
The CLI can still use native `minimap2`.

Facts found about other options:
- The `minimap2` Rust crate wraps C code (`minimap2-sys`). Its docs.rs build page for
  0.1.17 lists only `x86_64-unknown-linux-gnu`. No WASM build was found.
  Source: https://docs.rs/minimap2
- `kira-ls-aligner` is an early Rust aligner (v0.1.x). WASM status unknown.
- `yara-mapper` is built on C++ (SeqAn2). Ruled out.
- BWA is C. Earlier text said it is hard to build for WASM. That was wrong: https://github.com/pinbo/bwa-samtools-web
  (user-supplied link, read 2026-11-11) has BWA 0.7.17 and samtools 1.11 built with Emscripten, from biowasm
  patches. Only `bwa index`, `bwa mem`, `samtools sort` and `samtools index` are kept. `pthread` is removed
  (single-threaded). The author edited the C code so the last argument is the output file. 1 star,
  last push 2026-10-18, GPL-3.0. Not built or tested by us. `UNMEASURED`: its speed, memory and accuracy.
  `bwa mem` accepts paired files, so it could fill the paired-end gap of `rammap`'s streaming API.
  Weaknesses: separate C build (breaks the one-source goal), no known strong Nanopore mode, GPL-3.0
  licence effects on our project (the project is now GPL-3.0, see below; the legal effect of bundling
  GPL-3.0 code is not checked). Status: candidate second backend behind the `Aligner`
  trait. Not used in version 1. samtools is not needed in the core, because the caller reads the
  alignment stream and needs no sorted BAM.
- genomics-wasm (https://github.com/comparativechrono/genomics-wasm, user-supplied link, read
  2026-11-11 via the GitHub API and raw files): about 30 upstream tools compiled to WASM, created
  2026-10-10, last push 2026-10-24, 0 stars. Its own `LICENSE` is MIT (GitHub shows NOASSERTION).
  Relevant tools: minimap2 2.31, BWA 0.7.19, BWA-MEM2 2.2.1 (single-threaded), FreeBayes 1.3.6 and
  LoFreq 2.1.5 (variant calling), BCFtools and SAMtools 1.16 and 1.21 with plugins, mosdepth.
  `THIRD_PARTY.md` lists BWA as GPL-3.0 and FreeBayes dependencies with GPL-2.0-or-later code; the
  effect on our GPL-3.0 project was not analyzed. Stated limits: whole-genome performance is "not
  established"; FreeBayes heap up to 2 GiB with no promise that whole-genome runs fit; minimap2 works
  only for inputs that fit in browser memory. Not built or run by us; no speed or accuracy numbers.
  Possible use: an independent in-browser cross-check (minimap2 then BCFtools). Not used in version 1.
- biowasm (https://github.com/biowasm/biowasm, MIT, 280 stars, not archived; created 2019, last push
  2026-11-10, checked via the GitHub API on 2026-11-11). It is actively maintained. Its `tools` folder
  lists, among others: bcftools, samtools, htslib, minimap2, bwa, bowtie2, wgsim, fastp, seqtk, ssw,
  gffread and rammap. A `rammap` entry exists there. `UNVERIFIED`: what version it builds, how it is
  packaged, and whether it adds anything over our own `wasm-bindgen` build. A hosted biowasm `rammap`
  or `bcftools` could serve as a cross-check in the browser. Not used in version 1.
- MATPredict (`../MATPredict`) has no Rust aligner. It has C tools with WASM build
  scripts in `web/spike/native/` (miniprot, exonerate, losat, mmseqs2). `miniprot`
  aligns proteins to a genome. It does not map reads. Its own README says no Emscripten
  build was verified on that machine.

## 5. Variant calling

A custom haploid pileup caller in `qmutfind-core`.

1. Build a pileup from the CIGAR strings of the alignments.
2. Call a site when the alternate allele passes a minimum depth and a minimum read fraction.
   Both are user settings.
3. Call SNPs and indels.
4. Use separate default profiles for short reads and long reads.
   Nanopore reads have many indel errors, mostly in homopolymers.
5. Default thresholds are `UNMEASURED`. Set them from the benchmark (section 11).

Why custom: bcftools and FreeBayes are C and add WASM build problems. A haploid
caller is simpler than a diploid caller.

Classes:
- SNP.
- Small indel: 1 to 49 bp. The cutoff is a setting. Default 50 bp is a common
  structural-variant convention. It is not derived from this project's data.
- Large indel: 50 bp or more.
  - DECIDED (user, 2026-11-11): version 1 calls large indels from long-read alignments only.
    A split-read caller for short reads is deferred. Small indels are the main focus.
    A large indel in short reads usually appears as a split or clipped read, not as a
    CIGAR indel.
  - For short reads, the tool reports clipped or split-read evidence as
    "possible large indel, unresolved". It does not call it.
  - Label the large-indel class as low sensitivity in the output.

## 6. Effect classification

Needs a GFF3 with CDS features. UTR labels need UTR features.

Labels:
- SNP in CDS: synonymous, missense, nonsense (stop gained), stop lost, start lost.
- Indel in CDS: frameshift, in-frame insertion, in-frame deletion.
- Other: splice site, 5' UTR, 3' UTR, intron, intergenic.
- Splice site (changed 2026-11-11; the earlier text said "within 2 bp of an exon boundary"): the first 2
  and last 2 bases of an intron (the donor and acceptor positions). Exon bases next to a boundary keep
  their coding or UTR label. A "splice region" label for those exon bases is not in version 1.

Rules:
- An indel in a CDS is a frameshift when its length is not a multiple of 3.
  This includes 4 bp and 5 bp indels, not only 1 bp and 2 bp.
- Size class and effect are independent. A 1 bp deletion is a small indel and a frameshift.
- An indel outside a CDS has no frame. It gets a non-coding label.
- An indel that crosses an exon boundary or splice site is labeled `splice-disrupting`
  or `complex`. The tool does not guess a frame.
- UTR labels (DECIDED, user, 2026-11-11): use UTR features from the GFF when present. When absent,
  derive UTRs from the file: within each mRNA, the exon bases outside the CDS span are UTR.
  Bases on the 5' side of the start codon are 5' UTR. Bases on the 3' side of the stop codon are 3' UTR.
  Label derived UTRs `UTR (derived)`, so a user can tell them from UTRs given in the GFF.
- Limit, measured on Af293 (2026-11-11, `GCF_000002655.1` GFF3): of 9,823 mRNAs, 5,797 (59.0%) have a
  derivable UTR (3,563 both, 1,641 only 5', 593 only 3'). 4,026 (41.0%) have exons that end at the CDS
  ends, so no UTR can be derived. Median derived length: 268 bp (5') and 283 bp (3').
  A mutation just outside the CDS of such a gene is labeled intergenic, even if it lies in an
  unannotated UTR. Print this limit in the report.
- Introns (DECIDED, user suggestion, 2026-11-11): use intron features from the GFF when present. When
  absent, derive them as the gaps between consecutive exons of the same mRNA. Label derived introns
  `intron (derived)`. Splice sites are the first 2 and last 2 bases of each intron.
- Check on Af293 (2026-11-11): 21,586 derived intron instances (counted per mRNA, so introns shared by
  alternative transcripts are counted more than once) in 8,144 multi-exon mRNAs. Length min 1 bp, median 60 bp,
  max 17,227 bp. Ends: GT..AG 99.2%, GC..AG 0.6%, AT..AC 0.1%, other 0.1%. The 1 bp minimum is unexplained.
- If a GFF has no exon features and no UTR features, the tool labels the site `non-coding`. It does not claim a UTR.
- With no GFF at all, every effect label is `unannotated`. The tool still calls SNPs and
  indels and writes summary statistics. It cannot show gene context. It prints a warning
  at the start of the run.
- Later feature: project gene models from a related strain onto an unannotated assembly
  (for example with `miniprot`). Not in version 1.

## 7. Output

Main output: VCF.
- Fields: reference allele, alternate allele, depth, allele fraction, class, effect label.
- Each record has a `SAMPLE` field in `INFO`.
- DECIDED (user, 2026-11-11): the merged VCF has one record per sample per mutation, not a multi-sample
  column layout. A converter to the multi-sample layout may come later.
- `NSAMPLES` in `INFO` counts the samples that carry the same variant (same position and alleles).
  It makes no conclusion. A user can filter on it (for example `NSAMPLES=1`).
- Strand bias is a per-record field (see section 9).

Derived output: BED, written from the VCF.
The name column carries sample label, class and effect.
BED cannot hold alleles, so it is not the main file.

Per sample: one VCF. Composite: one merged VCF.
Comparison views between samples are deferred.

## 8. Report command

A separate command: `qmutfind report`.
- It reads the VCF or BED from an earlier run. It does not repeat mapping or calling.
- It draws locus views: gene model, mutation position, and effect.
- Version 1 selects mutations by genomic range, for example `chr3:100000-120000`.
- Design the selector to be extensible. Later selectors:
  - by feature (overlap with gene X in the GFF);
  - by computed property (for example "non-synonymous only").
- The report also shows the statistics from section 9.
- Version 1 draws locus views as self-contained SVG in an HTML page (gene models with derived UTRs and
  introns, mutations colored by effect). It also writes the VCF, BED and GFF for JBrowse 2 or IGV.
  Embedded JBrowse is a later option (decision 13 in section 12).

## 9. Statistics

Computed in `qmutfind-core`. Output: JSON (for programs), TSV (for reading), and the report.

Version 1:
1. Counts: total, by class, by effect label.
2. Ts/Tv ratio, with the transition and transversion counts shown next to it.
3. Substitution spectrum: 12 types and 6 strand-collapsed classes.
4. Indel length distribution, with the frameshift fraction.
5. Position: count and density (mutations per kb) per contig and in fixed windows.
   Default window 10 kb. This value has no data behind it. It is a setting.
6. Alternate-allele fraction histogram.
7. Callable genome and per-bp rate: the variant count divided by the bases that pass
   the depth and quality cutoffs. Print the cutoffs with the number.
8. Het-call rate: the fraction of calls with alternate fraction in a middle range
   (about 0.2 to 0.8, a heuristic). A flag, not an error count.
   Also flag sites with depth far above or below the median (about 2x and 0.5x, heuristics).
9. Observed versus expected non-synonymous to synonymous ratio, with a binomial
   confidence interval. Never print a bare dN/dS. With fewer than about 20 to 30 coding
   mutations the interval is very wide. A zero synonymous count gives an undefined ratio.
10. Share of indels in homopolymers and short tandem repeats, by repeat length.
11. Mutation clustering: distance to the neighbouring variant, with flags for groups
    within a set window.
12. Strand bias as a per-record VCF field. Not a hard filter at depth under about 20.
13. Trinucleotide context (96 channels), shown only above a minimum mutation count.
    DECIDED: the minimum is 100 mutations, adjustable (`UNMEASURED`, see open item 6).

Later: per-gene counts, genic versus intergenic enrichment, variants in repeats
(needs a repeat or mappability track, which the tool does not have).

Dropped: GC-content bias, CpG as its own statistic, RIP signature.
(RIP occurs in the sexual cycle. It is not expected in mitotic evolution.)

Source of items 7 to 13: an Opus subagent answering from general knowledge on
2026-11-11. It checked no sources. All thresholds above are heuristics.

## 10. Reference genomes and annotation

Input order:
1. Version 1: the user supplies FASTA and GFF3. Plain or compressed
   (`.gz`; `.zst` in the CLI).
2. Browser: pre-loaded genomes hosted on the same site as the app. This avoids CORS.
3. Later: fetch from NCBI by accession.
4. Later: FungiDB. UNVERIFIED: the user says it needs an API key. A key in browser code
   is visible to every user. So FungiDB fetch is CLI-only (key from an environment variable
   or config file), or each browser user types their own key.

Pre-loaded genomes:

| Strain | Notes |
|---|---|
| *A. fumigatus* Af293 | Standard reference. VERIFIED 2026-11-11 (NCBI Datasets API): `GCF_000002655.1`, ASM265v1, 29,384,958 bp, 19 contigs, RefSeq annotation. UNVERIFIED: that the GFF3 has CDS and UTR features. Check before building. |
| *A. fumigatus* A1163 | `GCA_000150145.1` (GenBank only, no RefSeq `GCF_` version), ASM15014v1, scaffold level, 29,205,420 bp, annotated (VERIFIED 2026-11-11). Also the user's reference for CEA10-background experiments. Web sources describe A1163 as a CEA10-lineage strain with an *A. niger* pyrG marker. One source calls it a CEA17 derivative. The conflict is unresolved (UNVERIFIED). |
| *S. cerevisiae* S288C | `GCF_000146045.2`, R64, 12,071,326 bp, 16 contigs, SGD R64-5-1 annotation (VERIFIED 2026-11-11). Small browser smoke test. ASSUMED strain. |

CEA10: NCBI has `GCA_051225625.1` (ASM5122562v1), a complete genome of 29,322,347 bp with no annotation
(VERIFIED 2026-11-11 for NCBI; EBI not checked).
Decision (user, 2026-11-11): map PRJNA1094460 to both references.
CEA10 assembly: no background differences, but no effect labels (tests the `unannotated` path).
A1163: annotated, so effect labels work, but with CEA10-versus-A1163 differences.
Comparing the two runs measures the background noise that the stand-in reference adds
(`UNMEASURED` until run).
A web result says Af293 is over 97% identical to A1163 (UNVERIFIED). Expect many background
differences between strain backgrounds.

CORS test (2026-11-11, one `HEAD` request with an `Origin` header):
- NCBI FTP genome FASTA (yeast): 200 with `Access-Control-Allow-Origin: *`. A browser can download it.
- NCBI Datasets API (405), FungiDB URL (404, a guessed URL), ENA portal API (403): the tests are inconclusive.
- Not tested: NCBI GFF files, streaming of large files in a browser.
- FungiDB (2026-11-11): the front page `https://fungidb.org/fungidb/` returns 200. A guessed downloads
  path (`/common/downloads/`) returns 404. I did not find the API-key requirement or the CORS
  behavior. Both remain UNVERIFIED.

## 11. Benchmark

"Matches bcftools" is not "correct". bcftools is also a caller. Two tests:

1. Simulated truth set.
   - Insert known SNPs, small indels and large indels into the Af293 genome.
   - Simulate Illumina-style and Nanopore-style reads. The simulator is not chosen yet.
     Shortlist from a web search (no source found that compares them for variant-caller truth sets):
     Illumina: wgsim, ART, DWGSIM, InSilicoSeq, Mason, NEAT (compared in Milhaven and Pfeifer,
     Heredity 2022, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9905089/; InSilicoSeq has the
     widest range of platform error models). Long reads: Badread (Nanopore), PBSIM3 (PacBio and Nanopore).
   - DECIDED (user, 2026-11-11): insert mutations with our own script (writes an exact truth VCF);
     simulate Illumina reads with ART and Nanopore reads with Badread. PBSIM3 and PacBio are out of version 1.
   - Availability checked 2026-11-11 (anaconda.org API): ART is in bioconda `linux-64`, latest 2016.06.05
     (an old release). Badread is in bioconda `noarch`, latest 0.4.2, and not in conda-forge.
   - Simulated reads are cleaner than real reads, so simulation overstates real accuracy. State this
     in the results. The real-data agreement tests cover the gap.
   - Built 2026-11-11 (test-first; `bench/make_truth.py`, 11 tests in `bench/tests/`, run with
     `pixi run test-bench`). The script reads FASTA (plain, `.gz`, `.zst`), designs non-overlapping
     variants (200 bp minimum spacing, 1,000 bp from contig ends, only on `ACGT` bases), and writes
     `<prefix>.mutated.fa` and `<prefix>.truth.vcf`. Same seed gives identical files.
     ASSUMED settings (no data behind them): small indel length 1 to 10 bp, large indel length
     50 to 500 bp, insertion and deletion equally likely, all variants at allele fraction 1.
   - Af293 truth set `af293_s1` (seed 1): 3,000 SNPs, 600 small indels, 100 large indels (counts are my
     choice). With 100 large indels, a sensitivity estimate for that class has a wide confidence
     interval. Files are in `/bigdata/stajichlab/jstajich/projects/qmutfind_data/bench/`.
     `bcftools norm --check-ref e` against the real reference reported no REF mismatch (3,700 records).
   - Finding: `bcftools norm` realigned 192 of the 700 indel records (27%), because random positions
     in repeats are not left-aligned. The truth VCF and every caller's VCF must be normalized with
     `bcftools norm` before they are compared.
   - Illumina reads: ART 2.5.8 (2016-06-06), HiSeq 2500 profile, 150 bp, 400 bp fragments (sd 50),
     paired, seed 1: 2,880,107 pairs, about 29x by my arithmetic (3:02 wall time, 46 MB). Files
     `af293_s1_art_1.fq.gz` and `af293_s1_art_2.fq.gz`. Paired files kept for the later pairing
     comparison (item 15).
   - FIRST MEASUREMENT (2026-11-11, `af293_s1` truth, ART 29x, `qmutfind call --preset short`, mates
     mapped as single-end, provisional thresholds `min_depth` 5 and `min_alt_fraction` 0.5 which are
     `UNMEASURED` and untuned; exploratory comparison script, both VCFs normalized with `bcftools norm`):
     SNP 2,985 of 3,000 found (99.50%), 7 false calls. Small indel 582 of 600 (97.00%), 12 false calls.
     Large indel 4 of 100 (4.00%), 0 false calls. Overall 3,571 of 3,700 found, 19 false calls.
     Run time 3:49 single-threaded, peak memory 476 MB, for 5,760,214 reads (release build).
     Limits: one seed, simulated clean reads, no `bcftools` comparison yet, thresholds not tuned,
     the comparison script is not yet part of the repository and has no tests. The low large-indel
     sensitivity is the design limit of short reads (decision: large indels from long reads only).
   - NANOPORE MEASUREMENT, first version (2026-11-11, `af293_s1_ont`: 102,791 reads, 822,777,529
     bases, about 27.7x by my arithmetic; `qmutfind call --preset long`; same provisional thresholds;
     tested comparison tool `bench/compare_calls.py`, both VCFs normalized with `bcftools norm`):
     SNP 2,999 of 3,000 found (99.97%), 10 false calls. Small indel 535 of 600 (89.17%), 42 false calls.
     Large indel 55 of 100 (55.00%), 1 false call. Overall 3,589 of 3,700, 53 false calls.
     Run time 3:34.6 single-threaded, peak memory 1.91 GB. Badread job time: 2,701 s for 4 parallel runs.
     Breakdown of the misses: large deletions 55 of 55 found, large insertions 0 of 45 found. Small
     insertions of 7 to 10 bp 59.3% found (70 of 118), small deletions of 7 to 10 bp 96.4% (108 of 112).
     Hypothesis (not proven): the caller counted insertions by exact inserted sequence, and per-read
     errors split the evidence across keys. A synthetic test reproduces that failure. The fix groups
     insertions by anchor and length and reports the per-position majority base. The rerun (v2) is
     the check on the real simulated reads.
   - NANOPORE v2 (2026-11-11, after the insertion-consensus fix; same reads, same thresholds):
     SNP 2,999 of 3,000 (99.97%), 10 false calls. Small indel 571 of 600 (95.17%), 42 false calls.
     Large indel 56 of 100 (56.00%), 1 false call. Overall 3,626 of 3,700, 53 false calls.
     Small insertions improved: 7 to 10 bp from 59.3% to 86.4% (102 of 118), 4 to 6 bp from 92.9% to
     96.9%. Large insertions are still missed: 1 of 45 found (the shortest). Large deletions 55 of 55.
     ART v2 is identical to v1 in every cell (no regression). Run times 3:35 (ONT, 1.97 GB peak) and
     3:13 (ART, 476 MB peak), measured while other jobs shared the CPUs.
   - NANOPORE v3 (2026-11-11; large insertions reported as symbolic `<INS>` with SVLEN, reads grouped
     by anchor within 20 bp and length within 10% (at least 5 bp); decision B of 2026-11-11; large
     indels compared with a tolerance of 50 bp in position and 0.8 in size ratio, both ASSUMED):
     SNP 2,999 of 3,000 (10 false), small indel 571 of 600 (42 false), large indel 100 of 100
     (1 false). Overall 3,670 of 3,700, 53 false calls. 45 symbolic records. At the two sites checked
     by hand: `NC_007194.1:28986` SVLEN 242 (truth 244), AO 26 of DP 29; `NC_007194.1:560432` SVLEN
     415 (truth 417), AO 23 of DP 25. Run time and memory of v3 were not recorded separately.
     ART v3 is identical to ART v2 in every cell and has 0 symbolic records (no short-read regression).
     The run times printed by the v3 job script are NOT v3 timings: a `sed` edit left the timing
     file path pointing at the v2 file. v3 timing and memory are unknown.
     CAUTION: the fix was designed after studying this same truth set (one seed), so 100 of 100 is a
     result on the data it was tuned on. It needs a second truth set with a new seed and new reads.
     The `bcftools` large-indel numbers below used the same tolerant rule.
   - BCFTOOLS BASELINE (2026-11-11; `minimap2` 2.31, `samtools` 1.25, `bcftools` 1.25; `bcftools call
     --ploidy 1 -mv`; both outputs normalized with `bcftools norm -f -m -any`; same truth set and the same
     tested comparison tool). Presets were not tuned. ART used `minimap2 -x sr` with the mates mapped as
     PAIRS and `mpileup -X illumina`. Nanopore used `minimap2 -x map-ont` and `mpileup -X ont-sup`
     (the plain `-X ont` profile contains `-I`, which skips indels; `ont-sup` is designed for very
     accurate reads, and the accuracy of the Badread reads was not checked). "Filtered" keeps
     `FMT/DP>=5 && FMT/AD[0:1]/FMT/DP>=0.5`, the same thresholds as `qmutfind`.
     ART, SNP: qmutfind 2,985 found / 7 false; bcftools raw 2,987 / 21; bcftools filtered 2,986 / 0.
     ART, small indel: qmutfind 582 / 12; bcftools raw 561 / 31; filtered 561 / 31.
     ART, large indel: qmutfind 4 / 0; bcftools raw 16 / 0; filtered 11 / 0.
     Nanopore, SNP: qmutfind (v2) 2,999 / 10; bcftools raw 2,999 / 28; filtered 2,999 / 0.
     Nanopore, small indel: qmutfind (v2) 571 / 42; bcftools raw 565 / 81; filtered 565 / 79.
     Nanopore, large indel: qmutfind (v2) 56 / 1; bcftools raw 67 / 0; filtered 66 / 0.
     Reading: SNP sensitivity is tied on both read types; bcftools with the filter makes no false SNP
     calls and qmutfind makes 7 (ART) and 10 (Nanopore). qmutfind found more small indels with fewer false
     calls on both read types. bcftools found more large indels than qmutfind v2 on Nanopore reads.
     Limits: simulated reads, one seed, `qmutfind` single-end vs `bcftools` paired on ART, untuned
     thresholds on both sides, run times not comparable (jobs shared 4 CPUs). `qmutfind` ignores base
     quality and mapping quality; this is a known gap that simulated reads may hide.
   - BWA-MEM2 + FREEBAYES BASELINE, ART reads (2026-11-11; `bwa-mem2` 2.2.1 AVX2, mates mapped as pairs,
     `samtools sort`, `freebayes` 1.3.6 `--ploidy 1` with default settings, `vcfallelicprimitives -k`
     then `bcftools norm -f -m -any`; same truth set and tested comparison tool). Raw: 18,342 records
     (17,747 SNP records, 14,769 of them false, mostly near-zero QUAL such as 3.8e-15 for 2 supporting
     reads), so raw FreeBayes output must be filtered. Filtered with `INFO/DP>=5 && INFO/AO/INFO/DP>=0.5`:
     SNP 2,978 of 3,000 found (99.27%), 0 false. Small indel 587 of 600 (97.83%), 3 false. Large indel
     2 of 100. Overall 3,567 of 3,700, 3 false calls. For comparison (filtered, same ART reads): qmutfind
     SNP 2,985 / 7 false, small indel 582 / 12; bcftools SNP 2,986 / 0, small indel 561 / 31.
     Timing (shared CPUs, not a speed comparison): `bwa-mem2 mem` 1:46 wall with 3 threads and 1.10 GB peak;
     whole map-and-sort stage 160 s; FreeBayes 2:20 wall. The FreeBayes peak memory printed by the
     job (21 MB) looks too low and was not checked, so it is not recorded as a result.
     `bwa-mem2 index` of the 29.75 Mb genome: 9.8 s, 715 MB peak, index files about 162 MB.
     Notes: BWA-MEM2 without `-R` gives the sample name `unknown` in the VCF. `vcfallelicprimitives`
     in vcflib 1.0.15 has `-k` (keep INFO) and no `-g`; it says `vcfwave` is preferred. My first run
     used `-g` and failed; the post-processing was redone.
   - CLI BACKENDS on the ART reads (2026-11-12; `qmutfind call` release build, mates joined into one file
     and mapped as single-end reads, `--threads 3`; same truth set, tested comparison tool, outputs
     normalized with `bcftools norm`; thresholds depth >= 5 and fraction >= 0.5 in the tool):
     `bwa-mem2` + built-in: SNP 2,986 of 3,000 (11 false), small indel 580 of 600 (11 false), large 0 of 100.
     `bwa-mem2` + FreeBayes: SNP 2,979 (0 false), small indel 580 (10 false), large 2 of 100.
     `rammap` + built-in (v3, same reads): SNP 2,985 (7 false), small indel 582 (12 false), large 4.
     Pairing matters for FreeBayes: the earlier PAIRED `bwa-mem2` + FreeBayes run had small indels 587 found
     and 3 false, against 580 and 10 false single-end. One seed; an observation only.
     `bcftools norm` realigned 592 of 3,571 FreeBayes records and 0 records from the built-in caller, so
     the FreeBayes output of the CLI is not left-aligned. Known gap: normalize it in the tool.
     Wall time (rough, 3 threads, one run each): `bwa-mem2` + built-in 1:52, `bwa-mem2` + FreeBayes 4:31.
     Peak memory printed by `time -v` (largest child process, not identified): 1.25 GB and 2.36 GB.
     Not built: `--reads2` for paired input; `.zst` reads are passed to `bwa-mem2` unconverted, which it
     cannot read (untested).
   - EVOL20 REAL-DATA COMPARISON (2026-11-12; Af293 RefSeq reference; parent SRR6277228, 3,037,384 pairs,
     mean depth 14.2x from `samtools depth`; EVOL20_3 SRR6277229, 1,427,487 pairs; four pipelines through the
     CLI, thresholds depth >= 5 and fraction >= 0.5; all outputs normalized with `bcftools norm`; mates:
     `rammap` as single-end, the others as pairs). Calls per sample (parent / EVOL20 / EVOL20-only):
     `rammap` + built-in 646 / 1,010 / 700; `bwa-mem2` + built-in 1,335 / 1,849 / 1,192;
     `bwa-mem2` + FreeBayes 874 / 1,315 / 923; `minimap2` + `bcftools` 88 / 79 / 36.
     The parent should be close to the reference, so most of the extra calls of the first three are
     probably false. This is an assumption: there is no ground truth for the parent.
     Diagnosis (parent, `rammap` + built-in): median DP at called sites 27 against a mean depth of 14.2
     (about 2x), and 387 of 646 calls have a fraction under 0.6. This fits reads from duplicated or
     repeated regions piling onto one locus. A true variant in a haploid should be near 1.0. The default
     `min_alt_fraction` 0.5 was provisional and `UNMEASURED`; the real data show it is too low for SNPs.
     Threshold sweep, parent calls at min fraction 0.5 / 0.7 / 0.8 / 0.9: `rammap` + built-in
     646 / 105 / 66 / 49; `bwa-mem2` + built-in 1,335 / 186 / 92 / 63; `bwa-mem2` + FreeBayes
     874 / 131 / 98 / 81. Overlap with the 88 `bcftools` parent calls falls as the threshold rises
     (`rammap` + built-in 76 / 56 / 50 / 43 of 88).
     Cost on the simulated truth set (found / false): ART SNPs 2,985 / 7 at 0.5, 2,975 / 1 at 0.7, 2,972 / 1
     at 0.8; ART small indels 582, 576, 551 at 0.5, 0.7, 0.8 (383 at 0.9); ONT SNPs 2,999 / 10, 2,991 / 0,
     2,933 / 0 at 0.5, 0.7, 0.8; ONT small indels 571, 431, 313. One fraction threshold hurts indels and
     long reads far more than SNPs, and the real-data excess is almost all SNPs (621 SNP and 25 indel
     calls in the parent).
     Candidate genes (RefSeq `NC_007198.1`): `ssk1` AFUA_5G08390 2,130,224-2,134,138 (-), AFUA_5G13710
     3,600,517-3,603,725 (-), `hrmA` AFUA_5G14900 3,851,274-3,852,983 (+). EVOL20-only variants in them,
     reported by all four pipelines: `NC_007198.1:2,132,142 T>TG` (ssk1) and `NC_007198.1:3,852,237 A>G`
     (hrmA). None in AFUA_5G13710 in any pipeline. The amino-acid effects are not known yet (no effect
     classifier). Only 6 EVOL20-only variants are reported by all four pipelines; 127 by three, 701 by two,
     1,044 by one. The gene names come partly from patent text and are unverified.
     Rough wall times (3 threads, one run each; parent / EVOL20): `rammap` + built-in 89 / 44 s;
     `bwa-mem2` + built-in 75 / 42 s; `bwa-mem2` + FreeBayes 168 / 86 s; `minimap2` + `bcftools` 134 / 84 s.
   - NEW DEFAULTS (2026-11-12; `CallerConfig::defaults_for`: depth 5, SNP fraction 0.8 for short reads and
     0.7 for long reads, indel fraction 0.5; accepted by the user; set from the sweeps above, so tuned on
     the same data they were measured on). Real data, calls old -> new (parent / EVOL20 / EVOL20-only):
     `rammap` + built-in 646 -> 75 / 1,010 -> 103 / 700 -> 61; `bwa-mem2` + built-in 1,335 -> 111 /
     1,849 -> 172 / 1,192 -> 114; `bwa-mem2` + FreeBayes 874 -> 127 / 1,315 -> 154 / 923 -> 100
     (`bcftools`: 88 / 79 / 36). Parent calls shared with the 88 `bcftools` calls: 53, 61 and 66. The
     `ssk1` insertion and the `hrmA` SNV are still EVOL20-only in all three pipelines.
     Simulated truth set, new defaults (found / false): ART `rammap` + built-in SNP 2,972 / 1, small indel 582 / 12,
     large 4; ART `bwa-mem2` + built-in SNP 2,973 / 0, small 580 / 11, large 0; ART `bwa-mem2` + FreeBayes
     SNP 2,979 / 0, small 580 / 10, large 2; ONT `rammap` + built-in SNP 2,991 / 0, small 571 / 42, large 100 / 1 false.
   - WASM ON REAL DATA (2026-11-12; Node 22.23.4, single thread, `pkg-node` built with SIMD128, module 566 KB;
     the Af293 gzip reference and the gzip reads streamed through `scripts/wasm-run.mjs`; the native CLI is the
     `rammap` + built-in path, release build, one run each): EVOL20_3, 2,854,974 reads: native 46.4 s wall,
     575 MB peak; WASM 59.4 s (1.5 s to index, 57.5 s to map), 696 MB resident. Parent, 6,074,768 reads: native
     89.0 s, 572 MB; WASM 123.8 s (1.3 s index, 121.9 s map), 654 MB. WASM is about 1.3 to 1.4 times slower
     than native. The WASM VCF is byte-identical to the CLI VCF on both samples (103 and 75 records).
     This is Node, not a browser (same V8 engine); no browser test yet (open item 14).
   - EFFECT ANNOTATION ON THE EVOL20 CALLS (2026-11-12; `qmutfind annotate` with the Af293 RefSeq GFF3;
     0.37 s and 55 MB for 103 records). EVOL20 labels: 79 intergenic, 6 synonymous, 5 3UTR_derived,
     4 frameshift, 3 missense, 3 intron_derived, 3 5UTR_derived. The two candidate-gene variants:
     `hrmA` AFUA_5G14900 `NC_007198.1:3,852,237 A>G` missense D304G (DP 5, AO 5, transcript XM_748074.2);
     `ssk1` AFUA_5G08390 `NC_007198.1:2,132,142 T>TG` frameshift (DP 6, AO 5). The `hrmA` missense agrees
     with the literature summary (itself partly from patent text, UNVERIFIED). Both are at the minimum depth,
     because EVOL20 has about 7x coverage. EVOL20 minus parent: 61 variants, 6 of them in genes: the two above,
     synonymous `AFUA_1G00480` T209T (`NC_007194.1:151,238`), synonymous `AFUA_8G04070` D179D and T204T
     (`NC_007201.1:884,282` and `884,357`, 75 bp apart), and an intron change in `AFUA_6G11620`. The other 55
     are intergenic and unexplained. `AFUA_5G13710` has no variant. Limits: UTRs and introns are derived,
     and 41% of Af293 mRNAs have no derivable UTR, so a variant in such a UTR is labeled intergenic.
   - Large insertions, direct check (2026-11-11, `minimap2` alignments of the same reads, not `rammap`):
     at `NC_007194.1:28986` (truth 244 bp) 26 of 33 primary reads show an insertion of 20 bp or more,
     with lengths 234 to 247. The most common length (242) is in 8 reads (24%). At `NC_007194.1:560432`
     (truth 417 bp) 24 of 29 reads show one, lengths 407 to 420 and one at 492; no length is shared
     by more than 5 reads. The caller groups reads by exact insertion length, so the support splits
     and stays under 0.5. Two sites only. Supports the length-spread explanation. Does not prove it
     for all 44 missed insertions.
   - Bugs found by running on real data, each fixed test-first: (1) the Af293 reference has IUPAC codes
     (for example `M`, `R`); the caller wrote them into the VCF `REF` field and `bcftools norm`
     rejected the file. The caller now skips any variant whose reference footprint is not A, C, G or T.
     (2) the insertion splitting described above. Synthetic tests missed both because they had no
     ambiguity codes and no per-read errors inside insertions.
   - Nanopore reads: Badread 0.4.2, fragment length mean 8000 and sd 6000 (my choice), default
     identity and error models (not inspected). Speed measured at 1x: 3:55 single-threaded (3,642
     reads, 29.4 Mb). The 28x set is four runs of 7x with seeds 1 to 4, in one background job.
     Result file: `af293_s1_ont.fastq.gz`.
   - For my caller and for `bcftools`, compute sensitivity and precision by class.
   - Use this test to set the default depth and allele-fraction thresholds, per read type.
2. Real data, agreement only.
   - Compare calls with `bcftools`. Report calls shared and calls unique to each caller.
   - There is no truth set, so label this "agreement", not "accuracy".
   - Where the paper publishes mutations (for example known *cyp51A* changes), use them as a partial check. `UNMEASURED`: whether the papers below publish such lists.

Real-data sets (Illumina, paired; catalog records only, nothing downloaded):

| Study | Runs | Description | Reference |
|---|---|---|---|
| PRJNA757188 | 3 | LH-EVOL strains evolved from an Af293 parent | Af293 |
| PRJNA417720 | 2 | Experimental evolution of hypoxia resistance, Af293 | Af293 |
| PRJNA1347423 | 12 | Directed evolution for triazole resistance, A1163 (ΔKU80), with controls | A1163 |
| PRJNA1094460 | 13 | Evolution on high copper. Wild-type lineages started from CEA10 conidia (per a web search summary of Handelman et al. 2024, mSphere) | A1163 or CEA10 (open, section 10) |

Literature check, 2026-11-11 (web search summaries, not full-text reads):
- PRJNA417720: strain EVOL20 evolved from Af293 (20 passages at 0.2% O2). Reported: three
  non-synonymous mutations versus Af293, including a missense change in Afu5g14900 (`hrmA`).
  Two other genes named: Afu5g08390 and Afu5g13710. Part of this comes from patent text, which
  is not peer-reviewed. Verify against the mBio paper (PMC5040115 and PMC7396965).
  This gives a small partial truth set: 3 expected coding mutations.
- PRJNA1094460: Handelman et al. 2024, "Evolution of the pathogenic mold Aspergillus fumigatus
  on high copper levels identifies novel resistance genes", https://pmc.ncbi.nlm.nih.gov/articles/PMC11332344.
  Reported mutations: Pma1 L424I, Cpa1 A37V, and a Gcs1 mutation (the change was not given in the summary).
  The search did not confirm the accession PRJNA1094460 itself. The PopSeq title matches.
- PRJNA757188 ("LH-EVOL"): no source found under that name. UNVERIFIED what the lines are.
- PRJNA1347423: the search failed with an API error. Not checked.
- Mutation lists are a partial check only. They list mutations the authors chose to report.

For PRJNA1094460, CEA10-versus-A1163 differences appear in every sample (if A1163 is the reference).
A variant in both the control and an evolved sample is background.
A variant only in the evolved sample is a candidate new mutation.

Source of the candidate list: a keyword search of the PopSeq catalog
(`/bigdata/stajichlab/jstajich/projects/PopSeq/catalog/runs/746128.Aspergillus_fumigatus.runs.csv.gz`,
4,639 runs, 135 studies, 44 keyword matches). I did not read the papers.
Other matches not selected: PRJNA680589 and PRJNA670081 (CEA10 ΔpksP and ΔKU80),
PRJNA334306 to PRJNA334314 (A1159 ΔKu70 T0 and T60), PRJNA9521 (single-end).

PopSeq holds only Illumina and BGI/DNBSEQ WGS. It has no long-read runs.
Long-read candidates from the ENA portal API (2026-11-11). The listing was cut at 40 rows for
Nanopore, so it is incomplete. These are natural isolates with no known mutations. They give agreement
data, timing and memory data, and not accuracy. Depth is my arithmetic (base_count / 29.4 Mb).

| Run | Study | Strain | Platform | Bases | Approx. depth |
|---|---|---|---|---|---|
| SRR28035789 | PRJNA1078818 | ATCC 46645 | Nanopore MinION | 0.95 Gb | 32x |
| SRR28035400 | PRJNA1078810 | 47-4 | Nanopore MinION | 1.16 Gb | 40x |
| SRR28036164 | PRJNA1078823 | B5233 | Nanopore MinION | 4.42 Gb | 150x |
| SRR39591699 | PRJNA1494181 | Af269 | Nanopore PromethION | 9.22 Gb | 314x |
| SRR31066031 | PRJNA1081109 | UD1 | PacBio Revio | 2.82 Gb | 96x |
| ERR8084627 | PRJEB50298 | not given | PacBio Sequel | 8.65 Gb | 294x |

`UNMEASURED`: how ATCC 46645 relates to CEA10. A PopSeq study (PRJNA670081) resequenced both.

Download rule (lab standard): use `prefetch` with `parallel-fastq-dump` as in the global
instructions. Verify read counts against ENA. Keep data on `/bigdata`.

## 12. Open items

1. DONE (2026-11-11): the `rammap` spike passed (section 4). Follow-up tests are items 14 to 17.
2. DONE for accessions (section 10). Af293 GFF3 checked: no UTR features (item 18). Still open: check the A1163 and S288C GFF3 files.
3. DECIDED (user, 2026-11-11): own mutation-insertion script, ART (Illumina), Badread (Nanopore). Section 11.
4. DONE in part: long-read candidates listed (section 11). The Nanopore list was cut at 40 rows.
5. Read the full text of the real-data papers. PRJNA417720 and PRJNA1094460 have partial mutation lists. PRJNA757188 and PRJNA1347423 were not resolved.
6. DECIDED (user, 2026-11-11): the trinucleotide (96-channel) statistic is shown only at 100 or more
   mutations. The value is adjustable and has no data behind it (`UNMEASURED`). Below the gate, the plain
   12-type and 6-class counts are always shown. Revisit with the benchmark, for example by subsampling calls.
7. Merged VCF layout: DECIDED (one record per sample per mutation). Large indels from long reads only: DECIDED.
8. DECIDED: run PRJNA1094460 against both the CEA10 assembly and A1163. Still open: the A1163 versus CEA10 relationship (sources conflict).
9. FungiDB: API-key requirement and CORS are still unverified. Test with a key, supplied by the user.
10. DECIDED (user, 2026-11-11): the Cargo workspace lives at the top level of this directory, with the
    three crates in subfolders. The user ran `git init`. `.gitignore` written. Large data stays on
    `/bigdata`, outside the repository. No remote yet. No commit made yet.
11. DECIDED (user, 2026-11-11): pin tools with a top-level `pixi.toml` (use `[workspace]`, not `[project]`).
    Pin: Rust toolchain and `wasm32-unknown-unknown` target, `wasm-pack` or `wasm-bindgen-cli`,
    `minimap2`, `bcftools`, `samtools`, and the chosen simulators. UNVERIFIED: which of these exist on
    conda-forge and bioconda for `linux-64`. Check each one when writing the file.
    Checked 2026-11-11 (anaconda.org API): bioconda `linux-64` has `art`, `bcftools` (up to 1.25),
    `samtools` (up to 1.25), `minimap2` (up to 2.31). Bioconda `noarch` has `badread` 0.4.2.
    Conda-forge has `wasm-bindgen-cli` (0.2.114 and later, 0.2.114 matches the `rammap` lockfile),
    `wasm-pack` (0.14.0 to 0.15.1) and `rust`. Not yet checked: the stable `rust` version string and the
    `wasm32-unknown-unknown` standard-library package for conda.
12. DECIDED (user, 2026-11-11): a pre-loaded genome bundle is one gzip FASTA, one gzip GFF3 and a manifest
    (strain, NCBI accession, assembly name, download date, file sizes, SHA-256). Browsers decompress gzip
    natively; the CLI also accepts `.zst`. Bundles are hosted as static files on the same origin as the app.
    A build script downloads from NCBI, checks the files and writes the bundle. A test checks checksums
    against the manifest. Bundles stay out of git; the script is in git. `UNMEASURED`: gzip size of each
    genome (my estimate for Af293 is about 9 MB, from a typical ratio of about 3.3).
13. DECIDED (user, 2026-11-11): version 1 draws locus views as self-contained SVG in an HTML page, with a
    gene-model track (including derived UTRs and introns) and a mutation track colored by effect. It also
    writes the VCF, BED and GFF so a user can open results in JBrowse 2 or IGV. Embedded JBrowse comes
    later ("we will try JBrowse later"). `UNVERIFIED`: JBrowse 2 embedding details; not tested.
14. Test the `rammap` WASM build in a real browser (Chrome and Firefox): memory, speed, Web Worker use.
15. DECIDED (user, 2026-11-11): version 1 maps each mate of a pair as an independent single-end read, in
    both the CLI and the browser, so the two builds give the same VCF. Pairing is a later feature, to be
    justified by the benchmark: compare single-end and paired calls on the simulated truth set and on the
    EVOL20 data (`UNMEASURED` until run). `rammap` offers a library API that could add pairing later.
16. DONE (2026-11-11): the `rammap` library API builds and runs for `wasm32-wasip1` with no
    `wasm-bindgen` entry points of our own, and its output is byte-identical to native. Probe:
    `spike/rammap/api_probe_main.rs`, runner `spike/rammap/run_wasi.cjs`, expected output
    `spike/rammap/api_probe_expected.txt`. A synthetic 200 kb reference with injected changes:
    150 bp read with 1 SNP gives `150M`, NM=1; 2 bp deletion gives `80M2D78M`; 3 bp insertion gives
    `70M3I80M`; a 3 kb Nanopore-like read with 5 mismatches, 1 D and 1 I gives NM=7 on both strands;
    a random read is unmapped. Notes: (a) `Aligner::from_seqs` takes in-memory sequences and `map_seq`
    returns structured `Mapping` values (coordinates, strand, MAPQ, CIGAR ops), so the core need not
    parse SAM. (b) `M` covers matches and mismatches, so the caller compares read and reference bases
    (or requests MD). (c) `map_pair` exists in the library. Its doc comments say the paired path loses
    per-mate detail (no `SA` tag, no `full_cigar`, no supplementary flag). Relevant to a later pairing
    decision (item 15). (d) Not tested: the `wasm32-unknown-unknown` target with our own `wasm-bindgen`
    wrapper (the `rammap` crate itself already builds there, see section 4).
17. Test the multi-threaded WASM build (nightly Rust, `SharedArrayBuffer`). It needs special HTTP
    headers on the host, so check hosting constraints first.
18. DECIDED (user, 2026-11-11): derive UTRs from `exon` minus `CDS` when the GFF has no UTR features.
    The Af293 RefSeq GFF3 has no UTR features (31,739 exon, 30,114 CDS, 10,085 gene, 9,823 mRNA).
    See section 6 for the label, the measured coverage (59.0% of mRNAs) and the limit.
19. Measure full-size runs in WASM (all 625k Nanopore reads, all Illumina reads of one run) to replace the
    extrapolation in section 4.

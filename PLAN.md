# qmutfind: plan

Licence (user decision, 2026-11-11): GPL-3.0. `LICENSE` holds the standard text (35,149 bytes,
byte-identical to the copy in pinbo/bwa-samtools-web). The Cargo manifests declare `GPL-3.0-only`.
ASSUMED: "only" and not "or later". Change `license` in `Cargo.toml` if you want "or later".
`rammap` declares `license = "MIT"` in its workspace `Cargo.toml` (read 2026-11-11). I did not read its
`LICENSE` file. The legal effect of combining MIT and GPL-3.0 code is not checked here.

Status (2026-11-12). Works and is tested (125 Rust tests, 25 Python tests, 2 Emscripten tests, 19 Node tests;
CI passed on GitHub on 2026-11-12 for the commit before the effect annotation; later commits have not run in CI):
- `qmutfind call` with selectable backends: `--aligner rammap|bwa-mem2`, `--caller builtin|freebayes`,
  `--reads2`, `--threads`, `--min-depth`, `--min-snp-fraction`, `--min-indel-fraction`, `--min-mapq`. The
  `rammap` path runs the same `Engine` as the browser build.
- `qmutfind annotate`: the effect of each variant on the canonical transcript of a gene, for any VCF
  (`EFFECT`, `GENE`, `TRANSCRIPT`, `AA`, `NTX`). DECIDED (user, 2026-11-12, question 43): the simple tags and the
  rule that each gene is judged on its canonical transcript (the longest CDS), with `NTX` giving the number of
  transcripts. The labels are plain words (for example `missense`), not Sequence Ontology names
  (`missense_variant`); the user was told this and agreed to the tags as defined.
- `qmutfind stats`: JSON and TSV statistics from a VCF. `qmutfind report`: one HTML page with an SVG of the
  gene models and the variants of a genomic range.
- WebAssembly: `Session` (reference bytes, read chunks, VCF out) and the functions `annotate`, `stats` and
  `report`, built with `pixi run build-wasm-node` or `build-wasm-web`; `pixi run test-wasm` checks that the WASM
  output equals the CLI output for all four. Emscripten 3.1.64 is installed by `pixi run setup-emsdk`.
- `README.md` documents the commands. `examples/browser` is a page and a worker. `pixi run test-browser` runs it in headless
  Firefox (Selenium, SpiderMonkey engine) on a small synthetic data set: the VCF, statistics and report equal the native CLI
  output byte for byte (plain and gzip input), and a truncated reads file shows an error (3 tests pass; a mutation check that
  hard-coded the sample name made 2 fail, 2026-11-12). Full-size real data, 2026-11-12 (headless Firefox, one worker, reads files via the page's file inputs): EVOL20_3
  (SRR6277229, 2,854,974 reads, 2 gzip files) finished in 63.9 s (66 s wall) and gave a VCF byte-identical to the native
  CLI with current defaults (49 calls; native 40.5 s). Peak resident memory of all Firefox processes: 1,484 MB (sampled
  every 2 s, so a short peak could be missed). A test found a real bug: Firefox's `DecompressionStream` fails on
  multi-member gzip (the `parallel-fastq-dump` output has 4 members), so the page now uses the vendored `fflate` (MIT);
  a test covers it. Not measured: Chrome, Safari, the larger parent sample (6.07M reads), Nanopore data in a browser.
  Note: the older `engine_*.vcf` files in `evol20/` were made before the mapq 20 default and are not comparable.
Not built: the statistics that need depth, the reference or per-strand counts (callable-genome rate,
trinucleotide context, strand bias, expected N:S), the "possible large indel, unresolved" report for short
reads, base-quality filtering, the genomics-wasm BWA-MEM2 and FreeBayes builds, a real-browser test, a
second simulated seed. The defaults for thresholds were set from data and are tuned on the data they were
measured on (see section 11).
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
     the FreeBayes output of the CLI was not left-aligned at that time. FIXED later: `left_align_calls` is now applied to the
     FreeBayes calls in the CLI (`qmutfind-cli/src/main.rs`). No new `bcftools norm` count was measured after the fix.
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
   - STATISTICS ON THE REAL CALLS (2026-11-12; `qmutfind stats` on the annotated `rammap` + built-in calls at
     the new defaults). Af293 parent: 75 calls (50 SNP, 25 small indel), Ts/Tv 1.0 (25 transitions, 25
     transversions), 15% of calls with a fraction from 0.2 to 0.8, 32 calls within 100 bp of another call.
     EVOL20_3: 103 calls (68 SNP, 35 small indel), Ts/Tv 4.67 (56 transitions, 12 transversions), SNP
     classes C>T 33 and T>C 23 of 68, 37% of calls with a fraction from 0.2 to 0.8, 41 clustered calls.
     Reading (an inference, not a result): the parent looks like noise (Ts/Tv near what random errors give);
     EVOL20 is richer in transitions, but 34% of its calls are indels and 37% are at mid-range fractions,
     so noise is still present.
   - LOCUS VIEWS OF THE TWO CANDIDATE GENES (2026-11-12): `qmutfind report` pages for `hrmA` and `ssk1` in
     EVOL20_3 (3.8 KB each; the `hrmA` page shows 3 transcripts and 4 CDS blocks in its window and the one
     variant with the tooltip `NC_007198.1:3852237 A>G missense D304G AFUA_5G14900`). Written to
     `qmutfind_data/reports/`, outside the repository. Not inspected in a browser.
   - MAPPING-QUALITY SWEEP (2026-11-12; `rammap` + built-in, new thresholds, `--min-mapq` q; both mates
     single-end; the min depth stays 5). Real data, parent / EVOL20 / EVOL20-only calls: q0 75 / 103 / 61;
     q10 69 / 58 / 17; q20 61 / 31 / 9; q30 54 / 18 / 5; q40 37 / 13 / 3. EVOL20 Ts/Tv: 4.67, 2.00, 1.67, 1.00,
     0.75. `hrmA` (3,852,237 A>G) is kept up to q20 and lost at q30 and q40 (its depth falls below 5);
     `ssk1` is kept at every q. Simulated ART (found / false): SNP 2,972 / 1 (q0), 2,980 / 0 (q10 and q20),
     2,978 / 0 (q30), 2,977 / 0 (q40); small indel 582 / 12 (q0), 580 / 10 (q10 to q30), 580 / 9 (q40).
     CORRECTION: the earlier reading that the EVOL20 Ts/Tv of 4.67 "looks more like real mutations" is not
     supported. Ts/Tv falls as the mapping-quality minimum rises, so the transition-rich calls were mostly
     low-quality mappings, consistent with divergent repeat copies (an inference, not tested).
     DECISION (taken by the assistant while the user was offline, following the approved plan to add quality
     filters and re-measure): the default minimum mapping quality is 20 for short reads and 0 for long reads.
     Short reads: measured on `rammap` only; `bwa-mem2` MAPQ values were not swept. Long reads: on the simulated
     Nanopore set (seed 1) `--min-mapq` 0, 10 and 20 give the same result (SNP 2,991 found / 0 false, small
     indel 571 / 42 at q0 and 571 / 41 at q10 and q20, large indel 100 / 1), so the filter is neutral there
     and no real long-read data exists to show a benefit. The long-read default stays 0.
   - SEED 2, AN INDEPENDENT TRUTH SET (2026-11-12; `af293_s2`: 3,000 SNPs, 600 small indels, 100 large indels,
     seed 2; ART 30x, 2,880,680 pairs, seed 2; Badread 15.8x (58,757 reads, 470,224,613 bases, four runs of 4x,
     seeds 11 to 14; lower coverage than seed 1); calls with the defaults in force at that time, short
     mapq 20, SNP 0.8 and 0.7, indel 0.5, nothing re-tuned). ART (found / false): SNP 2,970 / 3, small
     indel 586 / 6, large 2 of 100. This matches seed 1. ONT (found / false): SNP 2,974 / 1, small indel
     555 / 234, large indel 99 of 100. The large-insertion fix holds on new data. The Nanopore small
     indels do not: 234 false calls against 42 on seed 1 (precision 70%). Coverage and seed are
     confounded (15.8x against 27.7x), so the cause cannot be separated with these data.
     Profile of the 234 false small indels: 219 deletions and 15 insertions, median size 1 bp, median DP 9
     (true small indels: 15), median supporting reads 5 (true: 12), 204 of 234 at a fraction from 0.5 to
     0.6. This is the Nanopore 1 bp homopolymer deletion error at low local depth.
     Minimum supporting reads for indels, applied offline to the existing calls (small indel found / false;
     large indel found): seed 1: floor 0 571 / 42, 5 571 / 19, 6 571 / 16, 7 571 / 13, 8 571 / 13, 10 568 / 13
     (large 100 throughout). Seed 2: floor 0 555 / 234, 5 552 / 126, 6 549 / 72, 7 536 / 38, 8 509 / 20
     (large 97), 10 426 / 6 (large 86). SNPs unchanged (2,991 and 2,974).
     DECISION (assistant, user offline): default minimum supporting reads for indels 6 for long reads, 0
     for short reads (a floor of 6 would remove every call at 7x coverage; `hrmA` has 5 reads). Chosen as a
     balance: no loss on seed 1, 6 true indels lost on seed 2. This is tuned on both seeds, so it has no
     independent check yet. Floors of 7 and 8 remove more false calls and lose more true ones.
     CONFIRMED with the real binary (release build, `--preset long`, defaults): seed 2 SNP 2,974 / 1, small
     indel 549 / 72, large 99 of 100 / 0 false; seed 1 SNP 2,991 / 0, small indel 571 / 16, large 100 of 100 /
     0 false. These equal the offline numbers.
   - INDEPENDENT CHECK OF THE EFFECT LABELS AGAINST `bcftools csq` 1.25 (2026-11-12; `bcftools csq -l -f ref -g
     RefSeq GFF3` in local, one-variant-at-a-time mode, which indexed 9,895 transcripts, 30,114 CDS and 0 UTRs;
     `bench/compare_annotation.py`, 35 tests; compared on my canonical transcript, matching csq's empty
     transcript field on intron entries by gene, with gene Name aliases from the GFF). Variants: the 3,700
     simulated truth variants of seed 1, annotated by `qmutfind annotate` and by `csq`. First run found
     three real flaws in my classifier, each fixed test-first: (1) the first overlapping transcript was
     reported, so a variant in one gene's UTR and another gene's CDS got the weaker label (118 variants);
     (2) an indel anchored outside a transcript but affecting it was labeled intergenic (3); (3) an in-frame
     insertion or deletion that gains a stop, loses the stop or loses the start codon was labeled only
     in-frame (18; csq: stop_gained). A fourth flaw was found on real data and fixed: a deletion whose anchor
     is the last base of a gene but whose deleted bases are intergenic was `unclassified`; deletions are now
     judged by the deleted bases only. Final result on the 3,700: effect agrees 2,086, disagrees 32, both
     silent (intergenic) 1,055, csq silent where I label a gene 4 (2 non-coding RNA, 2 complex), my UTR
     labels where csq says intron 512 (csq has no UTR features), not comparable 11. All 32 disagreements are
     my `intron_derived` (30) or `5UTR_derived` (2) against csq `splice_region`, which csq extends to 8 bases
     into the intron while my `splice_site` is the first and last two bases (a difference of definition, not
     an error; no new label was added). Coding SNPs with an amino-acid change reported by both: 1,553, all
     identical in residue, codon number and new residue. Real calls: EVOL20 103 variants, 16 agree, 0 disagree,
     8 UTR-versus-intron, 9 of 9 amino-acid changes identical; parent 75 variants, 20 agree, 0 disagree, 9 of 9
     identical. Statistics cross-check against `bcftools stats`: SNP and indel counts, transitions,
     transversions and Ts/Tv are identical on the truth VCF (3,000 / 700 / 1,023 / 1,977 / 0.52) and on the
     EVOL20 calls (68 / 35 / 56 / 12 / 4.67).
     Limits: `csq` is a second implementation, not ground truth, and shares my reading of the GFF. The
     Af293 GFF has no mitochondrial contig and no `transl_table`, so the genetic-code handling is untested
     on real data. None of the 9,823 transcripts has a first-CDS phase other than 0, so my code, which
     ignores the GFF phase column, was not tested on partial CDS (it would be wrong on a GFF that has them).
     `csq` with `-p a` (haplotype mode) writes `BCSQ=@pos` pointers and looks silent for those records;
     local mode (`-l`) was used. SnpEff was not run.
   - INDEPENDENT CHECK AGAINST SnpEff 5.4c (2026-11-12; database built in 13 s from the same RefSeq FASTA and
     GFF3 with `-gff3 -noCheckCds -noCheckProtein`, in the pixi environment `validate`, which holds only SnpEff
     (1.8 GB); the same comparison script with `--snpeff`, 39 bench tests). SnpEff derives UTRs and introns
     itself, so it checks my derived labels. On the 3,700 truth variants: effect agrees 2,620, disagrees 7,
     both silent (intergenic) 1,055, SnpEff silent where I label a gene 6, not comparable 12. Amino-acid
     changes reported by both: 1,548, all identical. My `5UTR_derived` agrees with SnpEff's 5_prime_UTR_variant
     for 290 variants (14 more are refined by SnpEff to a premature start-codon gain), `3UTR_derived` with
     3_prime_UTR_variant for 206, `intron_derived` with intron_variant for 175. The 7 disagreements are large
     deletions at exon boundaries where SnpEff reports exon_loss, inframe_deletion or splice_region and I
     report splice_site or frameshift (a difference of priority). The 6 records where SnpEff is silent have
     lowercase letters in the REF allele (soft-masked reference), a SnpEff behaviour, not mine. Together with
     the `csq` check, the derived UTR and intron labels are supported by two independent derivations.
     Real calls against SnpEff: EVOL20 103 variants, 24 agree, 0 disagree, 79 intergenic in both, 9 of 9 amino-acid
     changes identical; parent 75 variants, 34 agree, 0 disagree, 41 intergenic in both, 9 of 9 identical.
     Limits: both tools read the same GFF; neither is ground truth.
   - WASM PARITY OF annotate, stats AND report ON REAL DATA (2026-11-12, Node 22.23.4; the full Af293 RefSeq
     GFF3 and gzip reference, the 103 real EVOL20_3 calls): `annotate` in WASM took 0.48 s (225 MB resident,
     plus 0.29 s to load the files) and its output is byte-identical to `qmutfind annotate`. The `stats` JSON and
     the `hrmA` `report` page are byte-identical to the CLI output too. A one-off script, not part of the tests.
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
9. FungiDB (2026-11-12): MEASURED with curl, no real key. The organism list (search `GenomeDataTypes`) is public and
   gives the FASTA and GFF URLs (Af293: 10,130 genes, release 71). The FASTA download (`/a/service/raw-files/...`) returns
   401 "Valid API Key required" without a key. The server reads the `Authorization` header; a fake `Auth-Key` header is
   ignored. The user reports that the key works; which form (`Bearer <key>` or the bare key) is not recorded, so
   `qmutfind fetch` tries `Bearer` first and the bare key second. `fetch` is tested against a local mock server and, on 2026-11-12, live.
   Live run (key file `--eupathdb-key-file`): the FASTA and GFF of Af293 downloaded (9.1 MB and 1.5 MB gzip, `gzip -t`
   OK). Two real-site findings fixed in the code: the GFF URL in the organism list (`/common/downloads/`) answers 404 and
   the file is under `/a/service/raw-files/`; and the key file was in `EUPATHDB_KEY=value` form. Contents of the FungiDB
   Af293 files (release 71): 9 sequences, 29,420,142 bases, named `Chr1_A_fumigatus_Af293` ... `Chr8_...` and
   `mito_A_fumigatus_Af293` (31,765 bp, so a mitochondrial contig IS available here, unlike the RefSeq GFF); 9,840
   `protein_coding_gene`, 9,840 mRNA, 30,134 CDS (the organism list said 10,130 genes, which probably adds the 289 ncRNA
   genes; not checked). Contig names differ from RefSeq (`NC_007194.1` ...), so results on one naming cannot be
   mixed with the other. One SNP inside a CDS was annotated by `qmutfind annotate` with this GFF (missense, E21K in
   Afu1g00100); the amino-acid result was not checked against another tool.
   CORS: no `Access-Control-Allow-Origin` header on the organism list, on the raw-files download or on a preflight
   (OPTIONS answers 401 on raw-files). From this, a browser page on another origin cannot read FungiDB responses, so a
   key-paste field in the WASM page would not work. NOT run in a real browser. Browser users upload the files instead.
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


## 11b. Mapq sweep for the `bwa-mem2` paths (2026-11-12)

`bench/run_mapq_sweep.sh`; `bwa-mem2` with the mates as pairs; defaults otherwise (depth 5, SNP 0.8, indel 0.5); q 0, 10,
20, 30, 40; `--caller builtin` and `freebayes`. VCFs are in `qmutfind_data/mapq_bwa/`; all normalized with `bcftools
norm -f -m -any` before comparing. For FreeBayes, q0 is sent as 1 (so q0 = q1).

Simulated ART seed 2 (truth set `af293_s2`; found / false): SNP found 2,970 to 2,983 for every q and both callers; false
SNPs 3 (built-in, every q) and 1 (FreeBayes, every q). Small indel 588 found at every q; false 6 (7 for FreeBayes q0).
Large indel 0 of 100 (short reads, by design). The mapq filter changes the simulated result by at most 0.2% of SNPs:
simulated reads have no repeat-driven mismapping, so this data cannot show the filter's value.

Real data (parent SRR6277228 / EVOL20_3 SRR6277229 / EVOL20-only calls; EVOL20-only Ts/Tv; `hrmA` and `ssk1` kept):
| q | built-in | FreeBayes |
|---|---|---|
| 0 | 111 / 172 / 114; 12.1; yes, yes | 127 / 154 / 100; 14.4; yes, yes |
| 10 | 109 / 140 / 86; 10.0 | 124 / 139 / 86; 14.3 |
| 20 | 108 / 132 / 78; 10.4 | 121 / 129 / 77; 16.3 |
| 30 | 99 / 109 / 56; 11.3 | 110 / 106 / 56; 16.0 |
| 40 | 98 / 98 / 45; 8.7; yes, yes | 106 / 96 / 47; 13.0; yes, yes |
`hrmA` (NC_007198.1:3,852,237 A>G) and `ssk1` (2,132,142 T>TG) are EVOL20-only at every q on both callers.

Comparison with `rammap` + built-in (same real data, q20): `rammap` gives 61 parent calls, 31 EVOL20 calls and 9
EVOL20-only calls, against 108 / 132 / 78 for `bwa-mem2` + built-in. The mapq filter removes far fewer calls with
`bwa-mem2` than with `rammap`. At q20 only 5 of the 9 `rammap` EVOL20-only calls are in the `bwa-mem2` + built-in set.
The 73 EVOL20-only calls that `bwa-mem2` + built-in makes and `rammap` does not (55 SNPs, 18 small indels) have median
depth 8 (range 5 to 25), median allele fraction 0.80, and only 10% have a fraction of 0.99 or more. `bwa-mem2` +
built-in and `bwa-mem2` + FreeBayes share 58 of their EVOL20-only calls, so the extra calls depend on the aligner, not on
the caller. They cluster: 21 of the 30 on NC_007195.1 lie within 1,854,035 to 1,861,813 (about 8 kb), and 15 of the 22 on
NC_007196.1 lie within 1,250,029 to 1,281,771.

Reading (inference, not tested): point mutations in an evolved line would not sit in dense clusters at fraction 0.8, so
these extra calls are probably mapping artifacts or a local structural difference, not independent mutations. No truth
exists for the real data, so this cannot be confirmed. It does mean `bwa-mem2` is not a drop-in improvement here, and
that "more calls" is not "more sensitive". The default mapq of 20 was chosen from the `rammap` data; for `bwa-mem2` it
helps less (114 to 78 EVOL20-only calls) and the evidence for a better value is weak. NOT done: look at the reads in
the clusters; paired versus single-end `rammap` (PLAN section 13, item 3).

## 12. Held-out validation (2026-11-12)

Pre-registered in `qmutfind_data/heldout/PREREGISTRATION.txt`; settings frozen at commit `bc69776`; no re-tuning
afterwards. Seeds 3, 4, 5 (standard) and 6, 7 (hard: `--min-spacing 30 --homopolymer-fraction 0.5
--cluster-fraction 0.25`). Each seed: 1,500 SNPs, 400 small indels, 100 large indels on NC_007196.1, NC_007198.1 and
NC_007199.1 (11.8 Mb). ART 5, 10, 20 and 30x (mates mapped as single-end); Badread 6, 12 and 24x. Defaults in force:
short reads depth 5, SNP 0.8, indel 0.5, mapq 20; long reads depth 5, SNP 0.7, indel 0.5, mapq 0, indel floor 6.
Raw rows: `bench/results/heldout_raw.tsv`; pooled table with 95% Wilson intervals: `bench/results/heldout_summary.tsv`
(made by `bench/summarize_heldout.py`). The pooled intervals treat variants as independent; the per-seed range is in
the table (`sens_min`, `sens_max`). Only 3 and 2 seeds, so the seed-to-seed spread is poorly known.
UNVERIFIED: the commit the release binary was built from was not recorded in the run log.

Sensitivity (found / truth) and precision, pooled over the seeds, standard mode:
| Reads | Class | Sensitivity | Precision |
|---|---|---|---|
| ART 5x | SNP | 0.542 (2,441/4,500) | 1.000 |
| ART 10x | SNP | 0.966 | 1.000 |
| ART 20x / 30x | SNP | 0.991 / 0.992 | 1.000 (0 false in 4,460 and 4,462 calls) |
| ART 5x | small indel | 0.515 | 0.984 |
| ART 10x | small indel | 0.932 | 0.987 |
| ART 20x / 30x | small indel | 0.975 / 0.975 | 0.986 / 0.985 |
| ONT 6x | SNP / small indel | 0.647 / 0.311 | 0.984 / 0.947 |
| ONT 12x | SNP / small indel | 0.977 / 0.820 | 0.999 / 0.920 |
| ONT 24x | SNP / small indel | 0.996 / 0.948 | 0.999 / 0.977 |
| ONT 12x / 24x | large indel | 0.920 / 1.000 | 1.000 |

Hard mode (2 seeds): ART SNP 0.948 at 10x, 0.991 at 20x, 0.993 at 30x, each with 1 false call; ART small indel 0.951,
0.981, 0.983 at 10, 20, 30x. ONT 24x: SNP 0.998 with 41 false calls (precision 0.986; standard mode had 3), small indel
0.935 (precision 0.932), large indel 0.990. ONT 12x hard: small indel 0.796, precision 0.874. The extra false SNPs in
hard mode were not diagnosed (the hard settings add clustering and homopolymer indels; the cause is not known).

Reading:
- At or above 20x, short reads give about 0.99 SNP and 0.975 small-indel sensitivity with precision 0.985 or higher.
  This agrees with the tuning-set results (PLAN section 11), so the thresholds did not overfit the first two seeds to a
  visible degree for ART.
- At 5x, sensitivity is about 0.54 for both classes. The minimum depth of 5 explains this: the median depth at 5x is
  near the threshold (inference, not tested by lowering the depth). Do not use the defaults below about 10x.
- Short reads find almost no large indels (0.3 to 2%), by design (PLAN decision: large indels from long reads only).
- Nanopore needs about 12x for SNPs (0.977) and 24x for small indels (0.948). At 6x, SNP sensitivity is 0.647.
  Nanopore small-indel precision is the weakest result (0.920 at 12x standard, 0.874 at 12x hard).
- Every number here is from simulated reads. Real reads have repeats and systematic errors that ART and Badread do not
  model, so real-data precision is expected to be lower. No real-data truth set exists yet.

## 12b. The site as three steps (2026-11-12)

Requested by the user: the region box is not useful before the variants are known. The page `site/index.html` now has
three steps that carry state: (1) Call (reference + reads -> VCF, kept in the page); (2) Summarize (statistics table and
a genome-wide plot; uses the VCF from step 1 or a VCF file the user chooses; takes contig lengths from the VCF `##contig`
lines or the reference); (3) Draw (annotate with the GFF, a table of variants with effect, gene, protein change, AF and
DP, a filter by effect, and a Draw button per variant that draws the locus with a window of 2.5 kb on each side).
Core `genome_plot_svg` (13 tests), `qmutfind plot`, and the WASM `genome_plot` give identical SVG (Node parity tests).
The browser tests (19) compare every output with the native CLI. Mutation checks caught: step 2 ignoring the carried
VCF, a wrong window, and (after adding a test) a stale annotated VCF surviving a new run.
The genome plot: top panel, calls per bin; bottom panel, the allele fraction of each call (Manhattan-style; there is no
p-value, so the y axis is the allele fraction, not -log10 p). Calls with no `AF` go on a separate lane.
Real-data check in Firefox (not an automated test): EVOL20_3 VCF (31 variants at the current defaults) with the RefSeq
genome and GFF: step 2 and step 3 ran, and Draw on the `hrmA` variant drew the missense change in AFUA_5G14900's
neighbourhood. NOT done: the window is a fixed size around the variant, not the gene extent; a whole-genome VCF with
thousands of calls was not tried in the browser; Chrome and Safari were not tried. `examples/browser` was removed:
`site/` is the only browser UI now.

## 12c. Genome menu, impact labels and where the browser can download from (2026-11-12)

- Impact labels: the locus drawing now prints the protein change (for example D304G) next to each variant, or a short
  effect name (`frameshift`, `intron`, `5UTR`) when there is none. A label that would cover another label or circle on
  its level is left out of the drawing. The page under the drawing lists every variant of the region in a table, so
  nothing is lost. The axis labels at both ends are anchored inside the plot.
- Genome menu (`site/genomes/genomes.json`, made by `bench/make_genome_list.py` from the NCBI RefSeq fungi table of
  2026-11-12): 155 genomes: all 134 RefSeq Chromosome or Complete Genome assemblies with annotation up to 70 Mb, the
  50 hand-chosen RefSeq genomes (16 Scaffold and 4 Contig level ones among them are kept by name), and two GenBank
  genomes (A1163 `GCA_000150145.1`, GFF3 present; CEA10 `GCA_051225625.1`, no GFF3 at NCBI, so `gff` is null and the page
  loads the reference only). 9 genomes are above 50 Mb (largest 67.2 Mb). Each has FASTA and GFF3 sizes from HEAD
  requests (all with `Access-Control-Allow-Origin`). The page downloads them in the browser and checks the sizes.
  Genomes above 50 Mb carry a warning: browser memory was measured only on a 29.75 Mb genome (about 0.7 GB in Node).
  Live check in Firefox: S. pombe (GCF_000002945.2) loaded in 1.1 s.
- FungiDB: MEASURED 2026-11-12 with a real key (key read from the key file, never printed). A GET with
  `Authorization: Bearer <key>` returns 200 but with no `Access-Control-Allow-Origin` header, and the preflight (OPTIONS)
  that a browser sends for an Authorization header returns 401 with no CORS headers. A web page therefore cannot download
  from FungiDB, with or without the key. The menu lists the 317 FungiDB organisms (FungiDB-71, public organism search,
  `site/genomes/fungidb.json` by `bench/make_fungidb_list.py`) and shows the `qmutfind fetch` command instead. A proxy
  would have to carry the user's key, which I did not build.
- Reads: the page downloads from ENA (HTTPS FASTQ, CORS `*`; 0.75 MB/s from the cluster). Alternatives tested: NCBI's
  on-the-fly FASTQ service (`trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq`) sends no allow-origin header; the AWS
  copy (`sra-pub-run-odp.s3.amazonaws.com`) is a `.sra` file (96 MB for SRR6277229, not FASTQ) with no CORS header. ENA is
  the only browser-readable source found. DDBJ was not tested.

## 12d. Gene names and products; readable plots in every theme (2026-11-12)

- The GFF3 `gene=` attribute (or, on a gene line, a `Name` that is not the locus tag) gives `GENENAME`; `product` on the
  mRNA, else on the CDS, gives `PRODUCT` (percent-encoded in the VCF so `;`, `=`, `,` and spaces cannot break the INFO
  field). The step 3 table has Name and Product columns and the locus drawing labels a gene `hrmA (AFUA_5G14900)`.
  On the RefSeq Af293 GFF 1,579 of 10,085 genes have a `gene=` name. The two candidate genes are named in the annotation:
  `hrmA` (AFUA_5G14900, product `protein hrmA`) and `ssk1` (AFUA_5G08390, product `response regulator`), which replaces the
  earlier note that these names came from unverified patent text.
- Plots: the browser test measures contrast in all five themes (text 4.5:1, marks 3:1, WCAG levels). Before the fix it
  found plot text at 2.3:1 in the quest theme (the plots used the page's foreground, white, on a grey panel), marks
  at 1.0 to 2.9:1 in several themes (the orange of the Okabe-Ito palette is 2.1:1 on white) and the locus drawing's
  missense colour at 2.25:1. Fix: one palette each for light, dark and mid panels, outlines in the panel's text
  colour, and the locus drawing (an iframe) receives its panel's colours and the matching palette after it loads and when the
  theme changes. Category colours are kept, but several are darker or lighter than the printed CLI report's.
  The standalone SVG and HTML files that the CLI writes keep the original palette.

## 12e. Web import and source lines (2026-11-12)

- Each input shows where its data comes from (local file, example, genome menu, ENA run, https address). Before this the
  inputs said "No file selected" while data from an example was in use.
- SRA run import: ENA's lookup service (`/ena/portal/api/filereport`, GET, `Access-Control-Allow-Origin: *`; HEAD returns 403)
  gives file addresses, sizes, MD5 sums, read count, layout, platform and organism; the page checks platform (Illumina or
  Nanopore only), the presence of FASTQ, the total size (limit 2 GB), the downloaded size and, after mapping, the read count
  (a paired run gives 2 x read pairs reads, because each mate is mapped separately). Live check in Firefox: SRR6277229 (paired,
  136.6 MB) and SRR28035789 (Nanopore, 920.7 MB, with a large-download warning) look up correctly; SRR0000000 is "not found";
  SRX1234567 is refused as a non-run accession.
- Risks considered. (1) Memory: a download is held whole in memory before mapping, so a large run can stop the tab; limited by
  the 2 GB cap and a warning above 500 MB; true streaming into the worker was not built. (2) Compression bombs: files loaded
  whole are capped at 1 GB unpacked; FASTQ is streamed and not capped. (3) Server access: only servers that send CORS headers
  work; the page says so. (4) Integrity: size and read count are checked; MD5 is not. (5) Privacy: the accession goes to ENA;
  no data goes to this site; private or access-controlled data is not supported. (6) Unsafe addresses: http, other schemes and
  addresses with credentials are refused. (7) Wrong data: platform and layout come from ENA's record, but a run of another
  organism than the reference is not detected. (8) Load on ENA: every user downloads from ENA directly; this site proxies
  nothing. NOT tested: Cancel during a real download; a download that stalls; Chrome and Safari.

## 13. Planned work (2026-11-12, agreed with the user)

1. Held-out validation (seeds 3 to 7): record the pooled results with confidence intervals. No re-tuning.
2. Mapq sweep for the `bwa-mem2` paths (`bench/run_mapq_sweep.sh`): real parent and EVOL20_3, simulated ART seed 2;
   `--caller builtin` and `freebayes`; q 0, 10, 20, 30, 40. For FreeBayes, q0 is sent as 1, so q0 equals q1.
3. Paired versus single-end mapping. Same data, `bwa-mem2` with the mates as a pair against the mates as two
   single-end files, then `rammap` single-end as now. Measure found/false on the simulated truth set and the call
   counts on the real parent. A `rammap` paired mode does not exist; whether it is worth building depends on this result.
4. Saved for later: validation against published mutation lists (PRJNA1094460, PRJNA1347423; the papers must be read).
5. Smaller items, one at a time, each with tests: base-quality filtering (off by default until measured), depth and
   strand statistics, a "possible large indel unresolved" report, Nanopore and the 6.07M-read parent in Firefox,
   `--caller bcftools` as an independent cross-check.
6. FungiDB files (contig names `Chr1_A_fumigatus_Af293`, with a mitochondrial contig): allowed as an input set. The
   source and release must be written next to every result (FungiDB release 71, fetched 2026-11-12 with `qmutfind fetch`).

### Redirects and Dropbox (2026-11-12)
A browser applies the CORS check to every hop of a redirect, and a page cannot read or change a redirect response. So
the page cannot follow a redirect whose first hop has no `access-control-allow-origin` header. A Dropbox share link
(`www.dropbox.com/s/...` or `/scl/fi/...`) is such a case: its 302 has no header (measured with `curl`). The same path on
`dl.dropboxusercontent.com` answers 206 with `access-control-allow-origin: *` and a working preflight (Range allowed).
The page therefore rewrites those links (`directAddress` in `site/app.js`; the `dl` and `raw` options are dropped). Other
hosts are not rewritten, and there is no general fix without a proxy, which was rejected. Not tested: a real browser
download of the 596 MB Dropbox file (headers only), and whether the `st=` part of a link expires.

### Read limit in the web import (2026-11-12)
A box "Keep the first this many reads per file" (default 2,000,000; 0 or empty = all) applies to SRA runs and to read
addresses. `downloadFirstReads` in `site/app.js` streams the file, unpacks gzip with fflate, counts four lines per read,
stops the download at the limit (`reader.cancel()`) and packs the kept reads again as gzip (level 3), so the page holds
only the kept reads. Both files of a pair are cut at the same read number. With a limit, a file above the 2 GB download
limit is allowed, because the download stops early. The expected read count after mapping is files x min(limit, reads
per file). The source line says "(first N reads of each file)".
- Not a random sample: the first reads of a run are used. Bias on real runs was not measured.
- A record is four lines, so multi-line FASTQ is not supported. Data after the cut point is never parsed.
- Live check (Firefox, headless, from the cluster): the 596 MB Dropbox R2 file with limit 100,000 loaded in 3.3 s and kept
  8.1 MB (about 160 MB expected at the default 2 M reads, extrapolated, not measured). Mapping a large kept set in the
  browser, and a limit applied to files chosen from disk, were not tested or built.

## 12f. Several strains: queue, catalog, combine and parent (2026-11-12; core, CLI, WebAssembly and page built)
Goal: process many strains (each a reads file R1 with an optional R2) one at a time in the page, name each strain,
keep a catalog, and combine the results, with an optional parent strain.

Built and tested:
- `qmutfind_core::combine` (`qmutfind-core/src/combine.rs`): input per strain = VCF text + optional coverage BED text.
  Rows = distinct (contig, pos, REF, ALT). Cells: call (AF:DP), `ref`, `low`, `unknown`, `other`. With a parent, status per
  strain: `new`, `inherited`, `lost`, `unknown`, `.`. A site where the parent has too little coverage is `unknown`, never
  `new`. `min_depth` (default 5 in the CLI, the caller's own minimum for short reads; NOT tuned) decides covered or low from the
  window mean depth. (It was 10 first; the review showed that this made `lost` impossible for a 5.9x strain.)
- CLI `qmutfind combine --strain NAME=VCF[,BED] ... [--parent NAME] [--min-depth N] --out-prefix P` writes
  `P.combined.tsv` and `P.summary.tsv`. WASM `combine_tables(names, vcfs, beds, parent, min_depth)` returns the two texts;
  a Node test checks that WASM and CLI give identical tables.

Not built (page): add-strain form (name, R1, optional R2, note, parent flag), queue (one strain at a time), catalog table
(name, files, reads, calls, time, settings, version), combined table view, export zip, import of an exported catalog,
IndexedDB persistence.

Decisions and reasons:
- The user adds strains one at a time (no file-name pairing). Settings are shared by all strains so the table is comparable.
- Memory: the engine drops each read after mapping (`Engine::map_records`); `Session::free()` ends the session. WASM
  linear memory does not shrink, so the page should stop and restart the worker after each strain (not built, not measured).
  Downloaded reads are held in `example.reads` until replaced; a batch must clear them after each strain.
- Tracking: an export zip is the durable record; IndexedDB is per browser and can be cleared.

Known limits and unknowns:
- Coverage windows are 10 kb in the page (`depth_bed(10000)`), so "covered" at a site is a window mean. A finer window may
  help; not measured.
- Indels from different strains must have the same left-aligned REF/ALT to share a row. CORRECTION (2026-11-12): before
  this date the built-in path did NOT left-align (only the FreeBayes path did); the review found it. `Caller::calls` now
  left-aligns every explicit indel, sorts again and keeps the call with more reads when two become equal (tests in
  `qmutfind-core/tests/sam.rs`). `<INS>` calls (anchor tolerance up to 20 bp), overlapping or nested indels are not merged.
- Review fixes in `combine` (2026-11-12): records with FILTER other than PASS or `.` are ignored; each ALT of a multi-allelic
  record is its own row; POS 0 and bad strain names are errors; a coverage contig that no VCF has is a warning.
- Parent status now uses an optional per-strain EVIDENCE VCF (`Caller::evidence_calls`: depth 2, fraction 0.1, two reads for an
  indel, no QUAL filter; `qmutfind call --evidence-out`, `Session.evidence_vcf()`). A cell with some reads for the ALT but no
  call is `weak`; a `weak` cell on either side makes the status `unknown`, never `new` or `lost`. Reason: on the real pair
  SRR6277228 / SRR6277229, 5 of 7 sites marked `new` had parent alternate reads (AF 0.18 to 0.40; review of 2026-11-12).
  The evidence thresholds (0.1, depth 2) are the reviewer's choice and were not tuned.
- Tests: `qmutfind-core/tests/combine_pipeline.rs` simulates a parent and a child through the caller with six sites of known
  truth (inherited, new, weak in the parent, lost, weak in the child, parent not covered). Simulated reads are clean; this
  checks the logic, not the thresholds.
- Per-strain time (about 80 s for 4 M reads) is an extrapolation from 58.8 s for 2.85 M reads, not a measurement.

### Parent status on the real pair after the review fixes (2026-11-12)
SRR6277228 (parent) and SRR6277229 (EVOL20_3), Af293, `rammap` + built-in, 1 kb coverage windows, evidence VCFs, `combine
--parent` with min depth 5, release build, both runs at the same time: parent 90.3 s and 585 MB, EVOL20 42.2 s and 634 MB.
- Before the fixes (10 kb windows, no evidence, min depth 10): 7 new, 22 inherited, 0 lost, 41 unknown.
- After: 3 new, 22 inherited, 6 lost, 39 unknown (evol20 has 31 calls, the parent 61).
- The five sites that the reviewer found with parent alternate reads (NC_007194.1:882403, NC_007195.1:1786423 and :1786426,
  NC_007195.1:2535794, NC_007196.1:4041535) are now `unknown`. The two clean sites (NC_007198.1:2132142 and :3852237) stay
  `new`. A third `new` site, NC_007196.1:2407906 TAG>T (EVOL20 AF 0.50, depth 6), was not in the reviewer's list and has weak
  support; it is not checked.
- The 39 `unknown`: 25 are parent calls where EVOL20 has weak evidence (likely shared, below EVOL20's thresholds at 5.9x), 8 are
  parent calls where the EVOL20 window is below the minimum depth, 6 are EVOL20 calls where the parent has weak evidence. A
  separate status for "parent call, child weak" may help; not decided.
- The 6 `lost` rest on a 1 kb window mean (at least 5) and no evidence record in EVOL20; site depth was not checked, so a site
  with 1 or 2 reads could be called `ref`. Treat `lost` as a lead, not as proof.
- The evidence files are large: 4.7 MB (0.32 MB gzip) for the parent and 16.7 MB (1.05 MB gzip, 180,449 records) for EVOL20.
  The low-coverage sample has more noise records at depth 2 and fraction 0.1. A floor of two alternate reads for SNP evidence
  may cut this; not tried.

### The strain section of the page (2026-11-12)
Built in `site/` (`index.html` section `#strain-section`, `app.js` "Several strains", `worker.js` `combine` step):
- Add a strain: name (required, unique, no control characters), reads R1 (required), R2 (optional), note. Names and notes are
  shown as text only (browser test with HTML in the name).
- "Run the waiting strains": one strain at a time, each in a NEW worker that is stopped when the strain ends (memory goes back
  to the browser; this also gives Cancel). A failed strain (bad reads) does not stop the queue. Cancel stops the running worker
  and marks the rest `cancelled`; they can run again. The 3-step buttons are disabled while the queue runs.
- Per strain the worker returns the VCF, a coverage BED with 1 kb windows, the evidence VCF, the qmutfind version, the read
  count and the time. Settings are shared: the reference and the read type of the page.
- Combine (needs two finished strains): parent select, minimum depth box (default 5), table with a filter (all, new in some strain,
  called in every strain, called in one strain), TSV download. Warnings go to the status line.
- Export: one zip with `catalog.json` (format `qmutfind-catalog`, version 1, reference name and size, settings, parent, one entry
  per strain with files, reads, calls, seconds, version), the strain files and, after a combine, `combined.tsv` and
  `summary.tsv`. Import reads such a zip back (size limit 500 MB for the zip and for the total unpacked size; at most 1000
  strains; every file named in the catalog must be in the zip; bad input changes nothing). Imported strains cannot run again.
- Not built: saving the catalog in the browser (IndexedDB). The evidence VCFs are large (see above); a quota and a gzip format
  would have to be measured first. Export is the record. Web sources (SRA, https) for strain reads are not wired to the list.

- LIVE CHECK of the strain queue on the real pair (2026-11-12, headless Firefox on the cluster, Af293 gzip reference, each
  strain with R1 and R2 as two reads files, parent first): the queue took 200.7 s in all (parent 6,074,768 reads, 61 calls,
  135.6 s; EVOL20 2,854,974 reads, 31 calls, 62.9 s; native runs of the same data took 90.3 s and 42.2 s). The combined table
  with the parent (70 variants) is byte-identical to `qmutfind combine` on the CLI VCF, BED and evidence files of the same
  data. Browser memory during the queue was not measured. Chrome and Safari were not tested.
- Test note: the browser tests read the state of the table in one script call, because the page rebuilds the rows and a held
  element goes stale (this made two tests flaky before the fix; three full runs passed after it).

## 12g. The page redesign: reference panel, sample set, results (2026-11-12)
User request: ask two questions first (reference genome; sample set), with the sources as a choice and rows added and removed
with + and -. The old step 1 (one reads input, Run) and the web-import panel were removed on the user's yes; a strain is a group
of one or more reads files, because the page does not use read pairs.
- Section 1, Reference genome: a dropdown (NCBI list, upload a FASTA and a GFF3, https addresses). One status line (`genome-status`)
  and the two source lines.
- Section 2, Sample set: read limit in THOUSANDS of reads per file (default 2000; 0 = all; a live echo says what it means), read type,
  and a table of sample rows. Each row: name, source dropdown (SRA run, web addresses, files from this computer), the reads fields of
  that source, a note and a minus button; `+ Add a sample` adds a row; one blank row exists at the start and blank rows are ignored.
  SRA rows have a Check button (ENA lookup and the old refusal messages); the lookup also runs when the sample's turn comes.
  Web downloads stop at the read limit (`downloadFirstReads`); the worker also applies the limit to every file (so local files
  are cut too, with no repacking). The run reports per file whether the limit cut it; the read-count warning is skipped then.
- Combine and catalog (the old strain section) and Results of one sample (a picker, plus a VCF/coverage file input for VCFs from
  other tools), then Summarize and Draw (steps 3 and 4). The example is one button that fills the reference and adds a sample.
- Coverage windows are 1 kb everywhere in the page now (they were 10 kb in the single-sample step).
- Tests: the browser tests were rewritten around the rows (122 pass). The step 1 tests now go through the helper `step1`, which
  fills one local sample 's1' and runs it. Mutation checks: no worker limit (local test fails), no duplicate-name check (fails).
- Not done: IndexedDB; per-sample read type (SRA sets it from the platform for its own run; the others use the page setting);
  checking the page at phone width and with a screen reader; Chrome and Safari.

## 12h. PRJNA1094460 (high-copper evolution) against A1163 (2026-11-12)
A real multi-sample test that is not EVOL20. Files in `bench/prjna1094460/` (sample table, download, calling job, analysis
script; the data are in `qmutfind_data/prjna1094460/`). This is not a unit-tested script; it is a record of one analysis.
- Design (ENA metadata): three genotypes, each with a control and evolved lines: CEA10 wild type (control `WT30`; lines A29,
  B23, C30, D30), `dCrpA` (control CrpA23; lines B23, C10, D15) and `dAceA` (control AceA6; lines B22, C22, D22).
- Data: the first 2,000,000 read pairs of each of the 13 runs (`fastq-dump -X`, a subsample of the first reads, not a random
  sample; each mate has exactly 2,000,000 reads), about 13.6x to 17.7x genome median depth. Reference: A1163
  `GCA_000150145.1` (55 scaffolds, 29.2 Mb, GenBank annotation with locus tags only).
- Run: built-in caller on all 13 in one SLURM job, 8 at a time: 5 min 13 s in all. Each sample has about 400 to 530 calls against
  A1163. `combine` of 13 samples with their evidence files: 6 s. The evidence files are 9 to 20 MB each.
- POSITIVE CONTROL (coverage, 1 kb windows): the two deletion strains lack the deleted gene. crpA: coverage is about 0 over
  DS499596.1:674-677 kb in all 4 crpA samples and normal in the other 9; the gene there is AFUB_036430, "copper resistance-associated
  P-type ATPase" (my first guess, AFUB_069550, was another copper ATPase and is covered normally). aceA: coverage is about 0 in
  DS499599.1:435-436 kb in all 4 aceA samples; the gene is AFUB_073740, "Copper fist DNA binding domain protein". Both mutants also
  lack DS499595.1:664-667 kb, AFUB_019720 "Ku family DNA helicase" (akuB/KU80): the mutants are in a KU80-deleted background.
- WT30 as the single parent does not fit the mutants (the 8 mutant samples share variants that WT30 lacks, including 4 SNPs near akuB).
  So each genotype was compared with its own control: new variants per lineage 14 to 21 for WT, 16 to 42 for crpA, 19 to 21 for aceA.
- THE PARENT ASSUMPTION IS NOT FULLY TRUE. Sites that are `new` in all lineages of a genotype cannot be independent events: all four WT
  lines share AFUB_054340 G>A (missense), AFUB_055510 T>C (missense) and AFUB_035300 (in-frame deletion) against WT30; all three crpA
  lines share the same C>T in AFUB_040930 (Ctr copper transporter family); all three aceA lines share the same frameshift in AFUB_099260.
  The simplest reading is that the control sample differs from the stock of the lines at these sites (or the lines came from one
  population that already had them). They are not evidence of adaptation until that is excluded. Variants private to one lineage
  are the cleaner candidates.
- AFUB_041460 (plasma membrane H+-ATPase Pma1) has coding changes in 8 lines of all three genotypes, at five different positions
  within DS499596.1:2,117,255-2,118,119. That is either a repeated target of selection or a mapping problem in this region. It was NOT
  checked (no alignment view, no second caller), and the published list of this study was not read, so no biological claim is made.
- Limits: no truth set; the calls were not checked with a second method; subsampled reads; the evidence thresholds and the 1 kb
  window were not tuned; the `lost` and `unknown` counts are large (9 to 51 and 125 to 213 per sample) and were not examined.

## 12i. Plot labels, gene arrows and a sortable variant table (2026-11-12)
- Contig labels of the variant plot and the coverage plot (`bands_and_labels` in `layout.rs`): bold, turned 90 degrees, right-justified
  (`text-anchor="end"` about the anchor) so that they hang BELOW `label_top` and never reach into the drawing. The coverage plot had a
  real overlap: its variant tick row uses y 254 to 262 and the labels were drawn at y 266; they now start at y 270. No two labels are
  closer than 12 px; when two collide the wider contig keeps its label (the others show nothing, but every band and label has a tooltip
  with the full name and length). A name over 24 characters is cut with an ellipsis. The picture is as tall as the longest label
  needs, and the legend of the variant plot follows the labels. Checked by eye on the 55-scaffold A1163 plots; tests in
  `tests/genome_plot.rs`, `tests/depth.rs` and `tests/common/mod.rs` (mutation checks: widest-first and the overlap rule).
- Locus drawing (`report.rs`): an arrow under each gene model (a line over the visible gene, an arrowhead at the end of transcription;
  at the plot edge, pointing out, when that end is outside the region) and chevrons along the gene line that point the same way.
- Variant table: every column header (not the Draw button) sorts: up, down, original order; numbers by value, text ignoring case with
  numbers in order, empty cells last in both directions, all variants sorted before the first 500 are shown, the sort stays when the effect
  filter changes and resets when new variants load.
- Build stamp: every file the page loads is named with a hash (`?v=`), because GitHub Pages lets a browser keep a file for 10 minutes and an
  old script could meet a new page. Not confirmed as the cause of the report that the reference dropdown did nothing; the deployed page
  worked in a fresh Firefox.

## 12j. All samples together (2026-11-12)
Request: summarize all variants together, to see which are shared, with the same filters as the combined table (called in one sample,
called in all samples, new, all variants), also in the annotation table.
- Core: `Combined::to_vcf` writes one VCF record per distinct variant: `NS` (number of samples with a call), `SAMPLES` (names, `%XX`
  escaped for `, ; = space %`), `NEW` (samples where it is new against the parent), mean `AF` and `DP`, `CLASS`, `##contig` lines
  and `##qmutfind_samples` / `##qmutfind_parent`. CLI: `qmutfind combine` also writes `PREFIX.combined.vcf`. WASM `combine_tables`
  returns it as the fourth text; the Node test compares it with the CLI file.
- Page: after a combine, the results picker has **All samples together**. The summary then shows the distinct variants, the rows
  "Samples combined", "Called in every sample", "Called in one sample only" and (with a parent) "New compared with the parent", and a
  table of calls, new, inherited, lost and unknown per sample. Annotation runs on the merged VCF; the variant table gets the columns
  Samples and Called in, and a "Shared by" filter (all, one sample only, every sample, new compared with the parent) that works with
  the effect filter and the sorting. A single sample has none of this. Editing or removing a sample makes the combine out of date and
  removes the choice.
- The AF of a merged record is the mean over the samples that call it; the coverage plot is not drawn for the merged view (the
  coverage of separate samples is not merged). Not done: colouring the genome plot by the number of samples.

## 12k. What a frameshift does to the protein (2026-11-12)
Request: a frameshift should say what truncation of the protein results.
- `annotate` (core `frameshift_protein`): the coding sequence is changed by the indel, extended with the exon bases of the 3' UTR (spliced,
  on the coding strand), and translated. The protein change is written like this: the first residue that differs from the reference, the
  residue it becomes, `fs`, and `*N` for the residues from the first changed one up to and including the new stop: `E3Dfs*4`.
  `E3Dfs*?` means no stop was found, not even in the 3' UTR. A first changed codon that is already a stop is `E3*fs`.
- New INFO tag `TRUNC=new/reference`: residues in the new protein over residues in the reference protein (without the stop), only when a
  stop is found. A value above 1 is an extension (the stop of the reference was lost in the shift).
- The page shows it next to the protein change: `E3Dfs*4 (protein 5 of 16 aa)`; the column still sorts by the change itself.
  The locus drawing uses the AA tag, so a frameshift is now labelled `E3Dfs*4` instead of `frameshift`.
- Tests: `tests/effects.rs` (fixtures c1, c2 and the new c6 and c7 with a stop that exists only in the 3' UTR, on both strands). Mutation
  checks: stop count, 3' UTR reading, minus-strand complement. Not handled: a frameshift near the start codon or across a splice site
  is judged on the spliced CDS only; an indel that spans a splice site is a `splice_site`, as before; only the canonical transcript.

## 12l. The Show dropdown of the Draw section (2026-11-12)
The Show dropdown (`#effect-filter`) is visible from the start (disabled, "Annotate first to filter the variants") and has two general
choices before the specific effects: All variants, and `Genic: in a gene (N)`: any variant that overlaps a gene (it has a GENE tag),
whatever its effect (missense, synonymous, frameshift, UTR, intron, splice site, non-coding exon, and so on). Intergenic and
unannotated variants are not genic. It works with the sorting and with the "Shared by" filter of the all-samples view.

## 12m. Which sample the variants belong to, the parent subtracted, and sorting in the combined table (2026-11-12)
- The chooser ("Choose what panels 3 and 4 show") and panels 3 and 4 are one coloured group. Each of the two panels has a "Showing" line:
  the sample, all samples together, or a VCF file, and the number of variants. What AF is: for one sample it is that sample's allele
  fraction (alternate reads divided by depth at the site); for all samples together it is the MEAN over the samples that call the
  variant, never a frequency of the whole data set. The note under the Show dropdown and the AF and DP header tooltips say the same.
  In the all-samples table, "Called in" lists each sample with its own AF (new INFO tag `AFS` of the merged VCF).
- "Subtract the parent" (a check box in the chooser): after a combine with a parent, panels 3 and 4 can show only the variants that the
  parent does not have: for one sample, the variants that are `new` in it; for all samples together, the variants that are new in any
  sample. The summary, the annotation table, the drawing and the VCF download follow it (the file is named `.new_vs_parent.vcf`).
  A variant for which the parent has weak evidence or too little coverage is `unknown` in the combine and is left out of this view,
  so the view is conservative. The parent itself has nothing to subtract; the check box is off when the combine is out of date.
- The combined table sorts by every header (up, down, original order; numbers by value; a sample column by the AF of its call; empty
  cells and the "no status" dot last in both directions); the sort stays when the filter changes and resets on a new combine.
- The locus drawing is titled with the sample (`Sample a: mutations in chrA:...`).

## 12n. Layouts, and sharing in the view of one sample (2026-11-12)
- Layout buttons next to Style: Standard, Wide (the page uses the whole window, with tables up to 75% of its height) and Phone (one column,
  bigger text and controls, full-width main buttons, form controls that fit the width, every sample as a card). `html[data-layout]`
  is chosen in the head before the page is drawn: `?layout=` first, then the saved choice (`localStorage` key `qmutfind-layout`), then
  automatic: phone when the window is 720 px wide or less, standard otherwise. An automatic layout follows the window; a click is a
  saved choice and stops that. Checked in Firefox only (headless; its smallest window is 500 px wide).
- One sample among several finished ones: its VCF gets the same sharing tags as the merged VCF (`NS`, `SAMPLES`, `AFS` and
  `##qmutfind_samples`), computed in the page from the finished samples. So the banner says "N variants, X only in this sample and Y
  shared with other samples", the summary has the rows Only in this sample / Shared with other samples / Shared by every sample, the
  table has the columns Samples ("3 of 13") and Called in (each sample with its own AF), the "Shared by" dropdown works (only in this
  sample, shared with at least one other sample, shared by every sample), and the annotated VCF download carries the tags. With one
  finished sample there is nothing to share and none of this is shown. A variant is the same in two samples when contig, position,
  REF and ALT are equal (the callers left-align indels, so the same event has the same form).
- New choice in the Shared by dropdown of the all-samples view: "Called in two or more samples".

## 12o. Pair-aware mapping with rammap (2026-11-12)
- Why: rammap maps each mate alone. A mate inside a repeat has mapping quality 0 and is dropped by the caller (min mapq 20), although its
  mate in unique sequence says where it belongs. `rammap` has `Aligner::map_pair`, but it returns PAF text without mate identity, so the
  pairing is our own, on top of single-read candidates.
- Core: `Alignment` has `score` and `primary`; `Aligner::align_candidates` gives primary and secondary alignments (`align` is unchanged);
  `pairing::choose_pair` (rule in the file header: proper pair = one contig, opposite strands, forward mate first, insert <= `max_insert`;
  best total score; a unique proper pair raises each mate to at least the quality of the other, capped at 60; no proper pair = old
  behaviour). `Engine::add_mate1/add_mate2/end_mates/set_max_insert`: chunks of the two files may split records anywhere; names must match
  (`/1`, `/2` and text after a space are ignored); different read counts are an error. `read_count` counts every mate.
- CLI: `call --paired` (needs `--reads2`; rammap path) and `--max-insert`. Default is off. The page and WASM do not use it yet.
- Tests: `qmutfind-core/tests/pairing.rs` (6), `pairing.rs` unit tests (7), two CLI tests. Mutation check: with the quality rescue switched
  off, the two paired-call tests fail.
- Evaluation, same binary with and without `--paired` (no truth set for the real data; scripts `bench/prjna1094460/compare_pairing.py`,
  `run_calls_paired.sh`, `run_calls_single.sh`, `bench/run_paired_heldout.sh`):
  - PRJNA1094460, 13 samples (2,000,000 pairs each, A1163): calls per sample 427 single, 429 paired (SNP 3069 -> 3115, indel 2487 -> 2466).
    Pairing gained 114 SNP and 25 indels and lost 68 SNP and 46 indels (pooled). Of the gained SNPs 88% are also called by another pipeline
    in the same sample; of the lost SNPs 82% are (so pairing also loses calls that others make). Recall proxy for background SNP 0.820 ->
    0.830, indel 0.928 -> 0.928 (the proxy is defined with the single-read set, so it favours it). Private calls per sample 26.6 -> 25.5.
    Jaccard with the bwa-mem2 built-in set: SNP 0.615 -> 0.630, indel 0.846 -> 0.849. bwa-mem2 still calls 1,376 SNP that paired rammap
    misses, so pairing closes about 4% of that difference (59 of 1,435).
  - Time and memory: 13 samples in 5 min 29 s paired against 5 min 13 s single; peak 600-665 MB in both.
  - Simulated ART pairs (seeds 3-5, 10x and 20x): SNP sensitivity 0.966 / 0.992 paired against 0.966 / 0.991 single; false positive SNP 0 in
    both; small-indel false positives 16 / 18 paired against 15 / 17 single. No regression, and no gain.
  - Reading: the effect is small. The cause of the larger gap to bwa-mem2 is not pairing; it was not investigated further here.

## 12p. Why bwa-mem2 + the built-in caller calls more SNPs than rammap (2026-11-12)
Scripts: `bench/prjna1094460/bwa_only_snps.py` (all 13 samples), `job_diag_wt30.sh`, `diag_reads.py` (WT30 only), and the diagnostic
`qmutfind-cli/examples/map_dump.rs` (one line per read: the rammap single-read alignment). Data in `qmutfind_data/prjna1094460/diag/`.
- The 1,376 SNPs that bwa-mem2 calls and paired rammap does not (13 samples): in 747 (54%) rammap has no alt evidence at all (fewer than 2 alt
  reads at mapq >= 20) although the bwa read depth is normal (median 10, AF 1.0); in 394 (29%) rammap has the alt at AF >= 0.8 but depth < 5;
  in 234 (17%) rammap has it at AF < 0.8. So the loss is mostly reads that rammap does not map or maps with a low quality, not the thresholds.
- WT30, 63 of these sites, 2,168 alt-carrying reads in the bwa-mem2 BAM, looked up in the rammap single-read output: 35% unmapped, 41% mapped at
  the same place with mapq < 20 (539 of them have mapq 11), 7% mapped elsewhere, 4% with an indel in the rammap alignment, 2% clipped, 10% mapped
  well (alt base counted). For the unmapped and mapq < 20 reads bwa-mem2 gives mapq >= 20 in about 90%, and their mate is not mapped by rammap
  within 1.5 kb either (2-4%), so mate rescue is not the reason and pairing cannot fix it.
- Two causes show in the examples. (1) Divergent reads: unmapped reads often have NM 9-24 over 150 bases against A1163 (CEA10-specific loci);
  bwa-mem2 aligns them locally (mapq up to 60), rammap's short-read preset finds no chain. At the 35 sites whose alt reads have median NM >= 6,
  32 lose at least half of the reads in rammap. (2) Repeat-like or multi-copy loci: some sites have 270-410 reads against a 15x window mean
  (the reference has fewer copies than the sample); bwa-mem2 piles them on one place at mapq 60, rammap gives mapq about 11 or leaves them
  unmapped. At the 28 sites with median NM < 6, 17 still lose at least half of their reads.
- Overall mapping is close: WT30 unmapped 3.80% (rammap, single reads) against 2.93% (bwa-mem2 primary reads not mapped); mapq < 20 2.48% against 1.91%.
  The difference sits in a small number of hard loci with many SNPs, which is where the extra bwa-mem2 calls are.
- What this does not show: whether the extra calls are right (no truth set). 55% of them recur in three or more samples, so most are
  CEA10-versus-A1163 background that the parent comparison removes. Sites seen in WT30 only (one sample, 63 sites) may not represent the other samples.
- Not tested: other rammap settings (seed length, minimizer window, chaining and local-alignment limits) on the lost reads, and a lower `min-mapq` together
  with them. A mapq 0 test earlier (WT30) did not close the gap (bwa-only calls 109 against 94).

## 12q. Saved results in the browser, and four review findings (2026-11-12)
Saved results (IndexedDB), `site/persist.js` and the "Keep my results in this browser" panel in the combine section:
- Saved after each finished sample, each combine, each removal of a sample and each import: the catalog zip (the same bytes as Export), the
  reference and gene files as loaded, and a small header (time, sample count, size, whether a combine existed). The reads are never saved.
  Keeping is on by default; the box is remembered in `localStorage` (`qmutfind-keep`). A reference or gene file above 300 MB is not saved.
- On the next visit the page only reports what it found ("Found saved results: ...") and shows Restore and Delete. Restore replaces the sample
  list (restored samples are "imported" and cannot run again), puts the reference and gene file back as if loaded from "this browser", and
  makes the combine again when one existed (it is deterministic). Saved results that were not restored are never overwritten: new results are not
  saved until the user restores, deletes, or presses "Save now (replaces the saved results)". An empty list never replaces saved results.
- Without IndexedDB (a private window, blocked site data, or `?nostore=1`) the page says so and everything else works. Write errors (quota) are
  shown and do not stop the page. The line shows the size of the saved data and `navigator.storage.estimate()` usage when the browser gives it.
- Not measured: quotas in other browsers; Chrome and Safari are untested. The data is per browser profile and the browser may clear it.
The four findings of the review:
- Gzip checksums: `site/gz.js` (`VerifyingGunzip`, with a CRC-32 in JS) checks the CRC and the length of every member, using fflate's `onmember`
  offset to read the trailer of the member before. Used for every gzip file that the worker reads and for the limited download. A file that is cut at
  the read limit is not read to its end and is not checked. A file whose trailing bytes are not a gzip trailer (padding after the last member) is
  an error now. Node tests: `tests/wasm/gz.test.mjs` (10; three of the checks fail when the comparison is removed).
- `FastqChunker` cap: an unfinished record may hold at most 64 MiB (`FastqChunker::with_limit`, default `DEFAULT_MAX_RECORD_BYTES`); past that it is an
  error that says the input may not be FASTQ. A 3 Mb record is fine. Browser test: 100 MB with no line break stops with "larger than".
- Read limit off the main thread: the unpacking, counting and repacking of `downloadFirstReads` now run in the worker (`fetchreads` step in
  `worker.js`, which needs no WebAssembly). Browser test: 300k reads from a web address; before the change the 50 ms heartbeat of the page was late by
  4,090 ms in total (longest gap 421 ms); the test now requires a longest gap under 400 ms and total lateness under 1.5 s.
- Total-size check: the reported bypass (look up with a limit, then set the limit to 0) was already closed, because the page looks the run up again
  at run time; a test now keeps it closed. A sample's files are now also counted together (a budget of the download limit for all its files; with a
  read limit, the kept packed bytes are counted). Browser test: three files that each fit and together do not.
- Tests now: 274 Rust, 42 Node, 183 browser (Firefox only). The bench tests were not run in this step.

## 12r. Fixes after the second review (2026-11-12)
A fresh review of 5dd1387..b0d139e (Opus, read-only) found one high and several medium problems. All are fixed except the two noted below.
- gz.js (high): `VerifyingGunzip` failed on good multi-member files when a chunk was larger than the 64 KiB tail (reproduced: 256 KiB chunks and
  a whole-file push failed). The tail is now cut to 64 KiB after the push, not before, so the whole chunk and the bytes before it are readable
  when fflate reports a member. An empty file now says "the gzip file is empty"; the checksum error mentions bytes after the last member.
  The reviewer's claim that exactly 64 KiB chunks fail when a member starts 0-9 bytes before a chunk edge could not be reproduced (tested
  member ends from -24 to +24 bytes of two chunk edges); the test stays as coverage. Node tests: 15 in `tests/wasm/gz.test.mjs`.
- Restore, import and combine while a sample runs are refused with "Wait until the samples have finished running" (before, a restore dropped the
  running sample and a combine took over the Cancel button).
- Two tabs: `saveNow` writes only if the stored header still has the `saved` time that this page last saw (checked inside the write transaction,
  `writeStore(entries, { expect })`). If another tab saved first, the page goes into the "Found saved results ... not saved" state. "Save now" still replaces.
- A restored reference or gene file now clears the file input of that slot, because `refFile()` prefers the input.
- `FastqChunker` scans only new bytes per push (kept scan position and newline count). 6 MB in 2 KiB chunks and 4 MB with no line break were 44 s and 90 s
  in a debug build before; both are far below 1 s now.
- Low: a write error inside a store transaction aborts it (all or none); the saved line no longer names a reference that was not saved (over
  the save limit); the save checks the unpacked catalog against the 500 MB limit that restore uses; the reference is written to the store only when
  it changed (counted in a test); `site/reads.js` (`keepFirstReads`, Node-tested, 9 tests) now holds the limited download: the download is bounded by the budget,
  the last flush is checked against `maxKept`, and a gzip file whose first chunk is one byte is still read as gzip.
- Test knobs added (like `download_limit_mb`): `?save_limit_mb=N`, `?catalog_limit_mb=N`.
- Not fixed: pairing can move a mate to a much worse secondary alignment and raise its mapq (pairing is off and stays unused until this is done: compare the pair
  score with the best unpaired score, as bwa does); CLI details of `--paired` (a value after it, `--max-insert` without it, memory when one file ends
  early). Not checked by the reviewer or by me: real browser chunk sizes, real quota pressure, Chrome and Safari.
- Tests now: 276 Rust, 56 Node, 190 browser (Firefox only).

## 12s. Checks before wider use (2026-11-13)
- Real gzip files through `site/gz.js` in Node (WT30 file `SRR28508724_1`, first 2 M lines, 207,545,380 bytes unpacked): plain gzip (1 member), bgzip (3,185
  members), two files concatenated (2 members) all unpack to the same bytes as `zcat` (md5 equal for the first two; the concatenated file has twice the bytes)
  and pass the CRC and length checks for chunk sizes of 16 KiB, 64 KiB, 1 MiB and 4 MiB; about 2 s per 200 MB. A file with 8 zero bytes after the last member is refused
  by fflate itself ("invalid gzip data"), as before gz.js existed; `gzip -t` accepts such a file with a warning.
- A real sample in headless Firefox, no read limit (SRR28508724, two files of 2,000,000 reads, about 360 MB gzip, reference A1163 plain 29.6 MB):
  4,000,000 reads mapped in 150.5 s and 152.5 s (two runs); peak resident memory of all Firefox processes 2.48 GB (sampled every 3 s, so a short peak can be missed).
  The result was saved 1.7 s after the run ("Storage used by this site: 15.9 MB"; the browser compresses it), the page offered it after a reload, and Restore took 0.6 s and
  gave the same VCF (469 lines). One sample per session only: memory across several samples was not measured.
- Not done: Chrome and Safari (no Chrome is installed here; a download of Chrome for Testing was not allowed in this session), a real phone, a real slow network.
- Speed against the commit before this work (5dd1387, built in a scratch copy; same machine, runs alternated; real sample SRR28508724, 4,000,000 reads, A1163):
  - WASM in Node (`scripts/wasm-run.mjs`, zlib unpacking, so this isolates WASM): old 128.85, 128.35, 128.74 s; new 127.94, 127.63, 127.97 s (mean 128.65 against 127.85 s,
    -0.6%); VCFs identical; peak memory 698-702 MB old, 701-710 MB new. The WASM file is 1,075,826 bytes now against 1,076,181 before.
  - Whole page in headless Firefox, no read limit: old 146.1, 146.5, 147.1 s; new 152.0, 150.0, 152.5 s (mean 146.6 against 151.5 s, +3.3%, +4.9 s); peak memory 2.47-2.50 GB for both.
    The extra time is the CRC-32 check of gzip files in JavaScript: fflate alone 1.31-1.43 s, with the check 1.84-2.03 s per 208 MB unpacked (+41% of unpack time).
  - CLI single-read calls, 13 samples: 5 min 13 s with the old binary and with the new one.

## 12t. Links out from the Gene column (2026-11-13)
- Page: "Link the Gene column to" (No links, FungiDB gene page, NCBI nucleotide), default No links, remembered in `localStorage` (`qmutfind-gene-links`). FungiDB:
  `https://fungidb.org/gene/<GENE>` (checked 2026-11-13: it answers with a redirect). NCBI: `https://www.ncbi.nlm.nih.gov/nuccore/<ACCESSION>`. Links open in a
  new tab with `rel="noopener noreferrer"`; the gene is always set as text and URL-encoded. No request goes to either site until a link is clicked.
- NCBI accession comes only from the GFF3: new INFO tag `NUCCORE` (core `annotate`, so the CLI, the WASM build and the page agree). Sources, in order: a `GenBank:` entry
  (any capitals) in `Dbxref` of the mRNA line (the Af293 RefSeq GFF has `GenBank:XM_001481640.2`), else a `transcript_id` of the mRNA line that has the shape of an
  accession (1-2 capitals, optional `_`, 5-12 digits, optional `.version`). Gene-line and CDS-line entries are not used. The A1163 GenBank GFF has neither, so
  with it the NCBI choice links nothing, and the page says so. FungiDB links use the gene (locus tag) as it is; a gene that FungiDB does not hold gives its not-found page.
- Tests: 7 in `qmutfind-core/tests/gene_names.rs` (two had passed without testing anything because their accessions were too short; fixed), 6 browser tests;
  mutation check: without the accession guard the "no accession" test fails. Tests now: 283 Rust, 56 Node, 196 browser.

## 12u. Sequence around the variant, comparison with snpEff, changelog (2026-11-13)
- `qmutfind-core/src/context.rs` (`sequence_context`, `context_text`; `SequenceContext` with `to_text` and `to_html`), CLI `context`, WASM `context`, and the page: Draw on a variant fills a monospace block under the locus
  drawing (setting: bases on each side, default 40). Rows: position, marker, reference, ref protein, variant, var protein. Exon bases are capitals; a minus-strand gene is reverse-complemented; an insertion adds a gap column; a
  deletion leaves gaps in the variant row; protein letters sit under the middle base of the codon; a frameshift reads into the 3' UTR; a variant that reaches the edge of the CDS shows a note instead of a protein. Tests: 11 core,
  4 CLI, 6 WASM-against-CLI (text and HTML), 5 browser (the page text equals the CLI text; alignment in a monospace block; flank setting; typed region hides the view; a gene named `G<b>1"x` is text). Seen on real Af293 data:
  a missense SNV (V29I), a 3 bp deletion, a 2 bp insertion with frameshift (A3Kfs*59) and a minus-strand missense (K164E) all read correctly.
- Comparison with snpEff 5.1 on 1,498,945 distinct variants (`bench/snpeff_validation/`, README "Validation of the annotation", report in `results/`): SNV agree 99.573%, convention 0.423%, differ 0.004% (60, all in tRNA genes);
  indel agree 96.812%, convention 0.743%, differ 2.445% (171, at exon edges and `complex`); protein change equal in 820,169 of 820,169 SNVs and 2,188 of 2,198 frameshifts; same transcript 99.950% (SNV). snpEff needed the
  `rna-`/`gene-` prefixes removed from `Parent` to link NCBI exons to transcripts (documented). qmutfind 14.5 s, snpEff 33.5 s for the set. Findings: `start_lost` is stricter than snpEff's (all 1,002 `M1L` are
  `initiator_codon_variant` in snpEff); 10 frameshifts near the gene end show a different residue and `?` for the reference residue (qmutfind display defect, not fixed).
- Earlier (2026-11-12) comparisons with `bcftools csq` 1.25 and snpEff 5.4c on 3,700 simulated truth variants were smaller; this run adds scale (all SNVs of 750 genes), protein-change checks and the gene/transcript check.
- `CHANGELOG.md` lists all features added so far. Tests now: 298 Rust, 62 Node, 201 browser, 74 bench.

## 12v. Confirming a knock-out or replacement allele: core and CLI (2026-11-13)
- Question from the user: introduced an hph replacement allele; where is it integrated and is it intact? Decisions: both read types with short reads first; CLI first, then a separate page; either construct form (cassette alone
  or a donor with arms), tested on simulated data and, if the marker can be found, on the public ΔcrpA and ΔaceA strains of PRJNA1094460 (the marker of those strains was not found: a web search did not identify it).
- Feasibility checked first: rammap gives two primary alignments for a read that crosses a genome-cassette junction, with exact coordinates (probe: 22 of 24,100 reads at 30x, junction at chr:50000 and chr:51500 as built).
- `qmutfind-core/src/confirm.rs`: `JunctionCollector` (hooked into `Engine::map_records`; ends are (contig, last base, side) so that a read of either strand gives the same junction; parts need 25 bases and mapq 1; neighbours on the read overlap by at most 50 or
  are separated by at most 50; clustering within 3 bases; at least 3 reads), `Confirmer`, `Report` (sites, cassette summaries, nine checks, verdict, text and JSON), `find_gene`, and a 24-mer comparison of the construct with the genome that warns about homology arms.
  Shared bases at a junction (microhomology) are reported and the flank that is not in doubt is given. CLI `confirm`; exit status 0/3/4.
- Tests: 11 scenario tests and 2 warning tests and 1 gene-lookup test in `qmutfind-core/tests/confirm.rs`, an opt-in sweep of 24 random worlds (passed: exact flanks, right orientation, confirmed, in all 24; 125 s), 5 CLI tests.
  Two first-run failures were useful: a one-base microhomology moved the reported flank (now reported as homology), and one test had a wrong coordinate after an earlier splice.
- Real reads (A1163, 4 M reads, dummy 2 kb construct, 95-100 s, 520 MB): ΔcrpA strain: `target_replaced` Pass (0% of genome depth inside AFUB_036430), no integration (correct); wild type: 105%, Fail. About 20 genome-to-genome junctions with 5 or more reads appear in each
  (repeats and artefacts); no check uses them, and the text report lists only the junctions that involve the cassette or lie within 5 kb of the target.
- Open: real transformant data; long reads; mate pairs as junction evidence; donor with arms (auto-strip by alignment of the construct to the genome); the page (`confirm.html`, worker step, WASM export); thresholds on real data.
- Tests now: 317 Rust (+1 opt-in), 62 Node, 201 browser, 74 bench.

## 12w. Allele confirmation page (2026-11-13)
- Core: `Confirmer::from_fasta` (genome FASTA and construct FASTA, with the name checks), `Report::junctions_tsv` and `depth_tsv`, and a fuller JSON (cassettes with gaps and variants, junction positions on the cassette, homology, target). The CLI uses the same functions.
- WASM: class `Confirmation` and `find_gene`; Node test: text, JSON and both tables equal the CLI, for two chunk sizes, with and without a target (4 tests).
- Page: `site/confirm.html` and `confirm.js` (same theme and layout scripts and the same worker as the main page; a `confirm` step in `worker.js`, which shares its read-feeding loop with `run`). Inputs: genome (file, https address, NCBI list), GFF3 (for a gene name), cassette (file or pasted
  sequence, with a name added when none), target (none, region or gene), reads (files and https addresses; a read limit uses the same worker download as the main page), thresholds. Result: verdict banner (words and a symbol; colour only as an addition), junction table,
  SVG depth plot (built with `createElementNS`, text set as text), checks table with chips, warnings, four downloads. The names from the files are always set as text (a test with a hostile contig name; a mutation to `innerHTML` makes it fail).
- Tests: 15 browser tests (inputs and links; a confirmed run whose downloaded report equals the CLI text; gene by name; wild type; low coverage; homology-arm warning; input errors; hostile names; cancel; readable in five styles; phone layout). A first draft of the check labels read badly in the banner ("Failed: Every part of the cassette is covered") and was renamed to noun phrases.
- Tests now: 317 Rust (+1 opt-in), 66 Node, 216 browser, 74 bench.

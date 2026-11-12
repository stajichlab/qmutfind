# qmutfind: plan

Status: design agreed in a question-and-answer review on 2026-11-11. No code exists yet.
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

Spike (do this before any other build work):
1. Build `rammap` for native and for `wasm32-unknown-unknown`.
2. Map one Illumina test file and one Nanopore test file to a small fungal genome.
3. Compare the output with native `minimap2` on the same input.
4. Record build outcome, time and peak memory.

Fallback if the spike fails: C `minimap2` built with Emscripten or WASI.
This keeps the browser build but loses the "pure Rust, one build" goal.
The CLI can still use native `minimap2`.

Facts found about other options:
- The `minimap2` Rust crate wraps C code (`minimap2-sys`). Its docs.rs build page for
  0.1.17 lists only `x86_64-unknown-linux-gnu`. No WASM build was found.
  Source: https://docs.rs/minimap2
- `kira-ls-aligner` is an early Rust aligner (v0.1.x). WASM status unknown.
- `yara-mapper` is built on C++ (SeqAn2). Ruled out.
- BWA is C and hard to build for WASM. Not used.
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
- Other: splice site (within 2 bp of an exon boundary), 5' UTR, 3' UTR, intron, intergenic.

Rules:
- An indel in a CDS is a frameshift when its length is not a multiple of 3.
  This includes 4 bp and 5 bp indels, not only 1 bp and 2 bp.
- Size class and effect are independent. A 1 bp deletion is a small indel and a frameshift.
- An indel outside a CDS has no frame. It gets a non-coding label.
- An indel that crosses an exon boundary or splice site is labeled `splice-disrupting`
  or `complex`. The tool does not guess a frame.
- With no UTR features, the tool labels the site `non-coding`. It does not claim a UTR.
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
- Locus views may use a static image first. A JBrowse view is a later option.
  The original plan named JBrowse. No decision was made on it.

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
    The minimum is `UNMEASURED` and needs a decision.

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

1. Run the `rammap` spike (section 4).
2. DONE for accessions (section 10). Still open: download each GFF3 and check that it has CDS and UTR features.
3. Choose the read simulator for the truth set (shortlist in section 11). Needs a user decision.
4. DONE in part: long-read candidates listed (section 11). The Nanopore list was cut at 40 rows.
5. Read the full text of the real-data papers. PRJNA417720 and PRJNA1094460 have partial mutation lists. PRJNA757188 and PRJNA1347423 were not resolved.
6. Set the minimum mutation count for the trinucleotide statistic.
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
12. DECIDED (user, 2026-11-11): a pre-loaded genome bundle is one gzip FASTA, one gzip GFF3 and a manifest
    (strain, NCBI accession, assembly name, download date, file sizes, SHA-256). Browsers decompress gzip
    natively; the CLI also accepts `.zst`. Bundles are hosted as static files on the same origin as the app.
    A build script downloads from NCBI, checks the files and writes the bundle. A test checks checksums
    against the manifest. Bundles stay out of git; the script is in git. `UNMEASURED`: gzip size of each
    genome (my estimate for Af293 is about 9 MB, from a typical ratio of about 3.3).
13. JBrowse integration: no decision yet.

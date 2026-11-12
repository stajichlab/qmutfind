# qmutfind

qmutfind finds mutations in one haploid fungal sample. It maps reads to a reference genome, calls
SNPs and indels, says what each variant does to the genes, and draws the result.
It is one Rust code base that runs as a command-line tool and as WebAssembly (WASM) in a browser.

Status: working and tested, not finished. Read [PLAN.md](PLAN.md) for the design, the measured results
and the open questions. Every number in this file comes from a measurement recorded there.

## Install

The tools are pinned with [pixi](https://pixi.sh). Run every command through `pixi run`.

```bash
pixi install
pixi run cargo build --release -p qmutfind-cli        # the binary is target/release/qmutfind
```

## Commands

`qmutfind call` maps reads and writes a VCF. One record per sample per mutation, with the sample name in INFO.

```bash
qmutfind call --reference ref.fa.gz --reads r1.fastq.gz --reads2 r2.fastq.gz \
              --sample mutant1 --preset short --out mutant1.vcf
```

| Option | Meaning |
|---|---|
| `--preset short\|long` | Illumina or Nanopore reads. Sets the aligner preset and the default thresholds. |
| `--reads2 FILE` | A second reads file. With `--aligner bwa-mem2` it is the second mate. With the default aligner it is mapped as more single-end reads. |
| `--aligner rammap\|bwa-mem2` | `rammap` (pure Rust, default) or `bwa-mem2` (short reads only, runs as an external program). |
| `--caller builtin\|freebayes` | The built-in haploid caller, or FreeBayes (needs `--aligner bwa-mem2`). |
| `--min-depth N`, `--min-snp-fraction X`, `--min-indel-fraction X` | Calling thresholds. Defaults: depth 5, SNP fraction 0.8 (short) or 0.7 (long), indel fraction 0.5. |
| `--min-indel-alt-reads N` | Call an indel only when at least N reads support it. Default 6 for long reads, 0 for short reads. |
| `--min-mapq N` | Ignore alignments with a lower mapping quality. Default 20 for short reads, 0 for long reads. |
| `--threads N`, `--tmpdir DIR` | For the external programs. The temporary directory is `--tmpdir`, then `$SCRATCH`, then `$TMPDIR`. |

Plain, `.gz` and `.zst` files are read directly.

`qmutfind annotate` adds the effect of each variant on the canonical transcript (the longest CDS) of a gene.
It works on any VCF, for example the output of `bcftools` or FreeBayes.

```bash
qmutfind annotate --vcf mutant1.vcf --reference ref.fa.gz --gff genes.gff3.gz --out mutant1.annot.vcf
```

It adds `EFFECT` (synonymous, missense, nonsense, start_lost, stop_lost, frameshift, inframe_insertion,
inframe_deletion, splice_site, 5UTR_derived, 3UTR_derived, intron_derived, non_coding, complex,
large_indel_in_cds, intergenic, unannotated), `GENE`, `TRANSCRIPT`, `AA` (for example `D304G`) and `NTX`.
UTRs and introns are derived from the exons and the CDS. Without `--gff` every record is `unannotated`.

`qmutfind stats --vcf FILE --out-prefix P` writes `P.stats.json` and `P.stats.tsv`: counts by class and effect,
Ts/Tv, the substitution spectrum, indel lengths, density per contig, an allele-fraction histogram, the share
of calls with a fraction from 0.2 to 0.8, and clustering.

`qmutfind report --vcf FILE --gff GFF --region CONTIG:START-END --out page.html` writes one HTML file with
the gene models and one circle per variant, coloured by effect.

`qmutfind plot --vcf FILE [--reference FA] --out genome.svg` draws every call along the contigs: calls per bin on top and
the allele fraction of each call below (a Manhattan-style plot). Contig lengths come from the VCF `##contig` lines, or
from `--reference`.

`qmutfind fetch --organism "Aspergillus fumigatus Af293" --out-dir genome` downloads the genome FASTA and GFF of one
FungiDB organism as `.gz` files. The organism list is public. The FASTA needs an API key from your FungiDB account
(Web Services Access). Put the key in the environment variable `EUPATHDB_KEY`, or in a file named by `--eupathdb-key-file`
(the file wins). A key file has one line per key: `KEY` for any site, or `HOST KEY` (for example `fungidb.org KEY`) for
one site. A line in shell form (`EUPATHDB_KEY=KEY` or `export EUPATHDB_KEY="KEY"`) also works. Blank lines and `#` lines are ignored. Run `chmod 600` on it; the command warns if others can read it. The key
is sent only to the `--base-url` host and never printed. `--base-url` changes the site (default `https://fungidb.org`). Tested against a local mock server and live
on Af293 (FASTA and GFF downloaded).

## WebAssembly

```bash
pixi run setup-emsdk        # Emscripten 3.1.64, only needed to build the genomics-wasm ports
pixi run build-wasm-node    # package for Node in pkg-node
pixi run build-wasm-web     # package for browsers in pkg-web
```

The package exports `Session` (the pipeline: reference bytes in, read chunks in, VCF text out) and the
functions `annotate`, `stats` and `report`. JavaScript decompresses the files before it passes the bytes.

```js
const session = new Session(referenceFastaBytes, 'short', 'mutant1');
for await (const chunk of fastqChunks) session.add_reads(chunk);   // a chunk can split a record
session.end_file();
const vcf = session.vcf();
```

Node: `node scripts/wasm-run.mjs --reference ref.fa.gz --reads r1.fq.gz --reads r2.fq.gz --sample S --kind short --out S.vcf`
prints the read count, the times and the memory in use.

Browser: the page in [site](site) has three steps. (1) Call: map the reads and call variants; the VCF stays in the page.
(2) Summarize: statistics and a genome-wide plot of where the variants are (`genome_plot`, the same code as
`qmutfind plot`). (3) Draw: add gene effects, list the variants, press Draw on one to see its locus. Each step uses the
result of the one before, or a VCF file you choose. `pixi run test-browser` runs the page in headless Firefox on a small
synthetic data set and compares every output with the native CLI. Chrome and Safari are not tested. Large real data sets
have been run only in Firefox (EVOL20_3, 2.9 million reads) in step 1.

## Tests

```bash
pixi run test          # Rust: the core and the CLI
pixi run test-bench    # Python: the benchmark scripts
pixi run test-wasm     # builds the WASM package, then checks it against the native CLI
pixi run test-emsdk    # the Emscripten install
pixi run test-browser  # the example page in headless Firefox (Selenium), compared with the native CLI
```

`test-wasm` checks that the WASM output is byte-identical to the native output for `call`, `annotate`,
`stats` and `report`.

## What has been measured

On a 29.75 Mb genome (*Aspergillus fumigatus* Af293), with simulated reads and a known truth set, and on
the real EVOL20 evolution experiment (PRJNA417720). See [PLAN.md](PLAN.md) section 11 for all tables.

- The WASM build gave the same VCF as the native CLI on two real samples (2.9 and 6.1 million reads), in
  1.3 to 1.4 times the time, with about 700 MB of memory (Node, one thread).
- On simulated reads, SNP sensitivity is 99% or higher for short and long reads. Small-indel sensitivity is
  about 97% for short reads and 95% for Nanopore reads. These simulated reads are cleaner than real reads.
- On real reads, repeats make false calls. The defaults were raised for that reason. On the real parent
  strain the calls fell from 646 to 75 with `rammap` and the built-in caller. `bcftools` found 88.
- Short reads cannot show large insertions in the alignment. Large indels are called from long reads only.

## Known limits

- Thresholds were tuned on the same data they were measured on. The parent-strain call count is a proxy, not truth.
- Nanopore small indels are the weakest class: on simulated reads at 15.8x, 72 of 621 small-indel calls are false even with the default floor of 6 supporting reads (16 of 587 at 27.7x).
- Mapping quality is filtered at 20 for short reads. That value comes from real short-read data only. Nothing was measured for long reads, so their default is 0. Base quality is not used.
- Statistics that need the read depth, the reference sequence or per-strand counts are not built.
- The `annotate` tags use plain words (`missense`), not Sequence Ontology names (`missense_variant`), and a gene is judged on its longest-CDS transcript only. Tools that expect Sequence Ontology terms or the `ANN` field will need a translation.
- Licence: GPL-3.0 (see [LICENSE](LICENSE)).

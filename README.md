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
| `--reads2 FILE` | A second reads file. With `--aligner bwa-mem2` it is the second mate. With the default aligner it is mapped as more single-end reads, unless `--paired` is given. |
| `--paired` | With `--reads2` and the default aligner: map the two files as pairs. A mate in a repeat takes the place that makes a proper pair with its mate (same contig, opposite strands, insert up to `--max-insert`), and the mate's mapping quality. Default off. |
| `--max-insert N` | The largest insert for `--paired`, in bases (default 1000). |
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
`NUCCORE` is the NCBI nucleotide accession of the transcript, when the mRNA line of the GFF3 gives one (a `GenBank:` entry in `Dbxref`, as RefSeq files do, or a `transcript_id` that is an accession); a GFF3 without either has no `NUCCORE` tag. In the web page, "Link the Gene column to" makes the gene a link to FungiDB (`https://fungidb.org/gene/GENE`) or to `https://www.ncbi.nlm.nih.gov/nuccore/ACCESSION`; a gene without an accession is not linked in the NCBI mode.

`qmutfind context` shows the sequence around one variant as aligned text rows (it is the text under the locus drawing in the web page):

```bash
qmutfind context --reference ref.fa.gz --gff genes.gff3.gz --variant NC_007194.1:300:G:A [--flank 40] [--html] [--out FILE]
```

```
NC_007194.1:300 G>A  gene AFUA_1G00100 (rna-XM_001481640.2, plus strand)  effect missense V29I
position    260       270       280       290       300       310       320       330
marker                                              ^
reference   GTTTTCAGCGTTCAAGGAGGTCCCGTGTATGCTGTTCTGCGTAGCGATGTTCTTTGGTTATATTGGATTCTTCAACCCCAT
ref protein   F  S  A  F  K  E  V  P  C  M  L  F  C  V  A  M  F  F  G  Y  I  G  F  F  N  P  I
variant     GTTTTCAGCGTTCAAGGAGGTCCCGTGTATGCTGTTCTGCATAGCGATGTTCTTTGGTTATATTGGATTCTTCAACCCCAT
var protein   F  S  A  F  K  E  V  P  C  M  L  F  C  I  A  M  F  F  G  Y  I  G  F  F  N  P  I
```

Exon bases are capitals and all other bases are small letters. `^` marks the changed column. The protein letter sits under the middle
base of its codon (a codon that an intron splits is translated across the junction). An insertion adds a column with `-` in the reference
row, and a deletion puts `-` in the variant row. After a frameshift the variant protein reads on in the new frame, into the 3' UTR if needed.
A gene on the minus strand is shown as it reads (the reverse complement of the reference). Outside genes the rows show the forward strand and
no protein. `--flank` is the number of bases on each side (default 40, at least 1). `--variant` is `CONTIG:POS:REF:ALT` in VCF form (an indel has
one anchor base). The command stops with a reason when the contig is missing or when REF does not match the reference. `--html` writes the same
text with the changed columns in `<span>` tags. The WASM function `context` gives the same text.

### Confirming a knock-out or a replacement allele

`qmutfind confirm` answers one question: *I introduced a marker cassette (for example `hph`); where is it integrated, and is it intact?*

```bash
qmutfind confirm --reference genome.fa.gz --construct hph.fa --reads r1.fastq.gz --reads2 r2.fastq.gz --preset short \
  --gff genes.gff3.gz --gene AFUA_1G00100 --out-prefix mutant1
```

`--construct` is a FASTA file with the sequence you introduced, **the cassette only, without homology arms**. `--target CONTIG:START-END` (1-based) or
`--gene NAME --gff GFF3` says which region the cassette should replace; without either, the tool still reports where the cassette is. Reads are
mapped to the genome plus the cassette as one more contig. The report is printed and written to `PREFIX.confirm.txt`; the same result is in `PREFIX.confirm.json`,
all junctions are in `PREFIX.junctions.tsv`, and the depth at each cassette base is in `PREFIX.cassette_depth.tsv`. The exit status is 0 when confirmed, 3 when not
confirmed, 4 when inconclusive (1 for an error); `--always-exit-zero` gives 0 after a report. `--min-junction-reads` (default 3), `--min-depth` (default 10) and
`--tolerance` (bases between a flank and the target edge, default 20) change the thresholds.

How it decides. A read that starts in the genome and goes on into the cassette aligns in two parts. The tool groups such split reads into junctions (the last
reference base next to the breakpoint and the side the read lies on; a read from either strand gives the same junction). Where the two parts share bases, the breakpoint is
uncertain by that many bases and the report gives the flank that is not in doubt. Each check is Pass, Fail or Inconclusive; a lack of coverage is never a pass.

| Check | Pass when |
|---|---|
| `integration_found` | at least `--min-junction-reads` reads join the genome to the cassette (Inconclusive when the genome depth is below `--min-depth` and nothing was found) |
| `both_sides_joined` | the main site has a junction on the left and on the right |
| `single_site` | there is one site; several sites mean extra insertions |
| `at_target` | both flanks are within the tolerance of the target edges (needs a target) |
| `target_replaced` | the inside of the target has at most 10% of the genome depth (no wild-type copy; needs a target) |
| `cassette_coverage` | every cassette base has at least `max(3, 0.2 x genome depth)` reads |
| `cassette_sequence` | no difference from the given cassette sequence is called (SNP or small indel) |
| `cassette_structure` | no read joins one part of the cassette to another (an internal deletion, a duplication or a concatemer) |
| `copy_number` | cassette depth is 0.6 to 1.5 times the genome depth (one copy per genome, for a haploid) |

The verdict is *confirmed* when nothing failed and nothing was inconclusive. The report also lists the orientation of the cassette (forward when its start is joined to the left flank),
the number of genome bases that it replaced, and warnings. A warning says that the construct shares sequence with the genome (24-base words found in both): usually these are homology arms,
and reads from them map to the genome and to the cassette equally well, so remove the arms from the construct file.

**What was tested.** Simulated 30x reads from a 120 kb genome with a 2 kb cassette: an intact replacement (24 random genomes, both orientations, replaced lengths 300 to 2,800 bases: the flanks were
exact and the verdict *confirmed* every time), a reversed cassette, a point change and a deletion inside the cassette, a missing stretch inside the cassette, a second copy elsewhere, a half wild-type
population, an insertion away from the target, no insertion, and 2x coverage (inconclusive). On real reads of *A. fumigatus* A1163 (4,000,000 reads, about 95 s and 520 MB) with a dummy construct that is not in the reads:
the Δ*crpA* strain gave no integration and 0% coverage inside AFUB_036430 (`target_replaced` passed), and the wild type gave 105% (it failed), as expected.
**Not tested:** a real transformant with its real marker; long reads (`--preset long` runs, but no thresholds were checked on long reads); reads whose two mates span a junction (mates are mapped as single reads and
only split reads count); diploid or heterokaryotic strains beyond the half wild-type simulation; constructs with homology arms (the tool warns but does not remove them). The thresholds were not tuned on real data.
For short reads the tool shows the junctions and the coverage; it cannot show in one piece that the whole insertion is contiguous, which only long reads can.

`qmutfind stats --vcf FILE --out-prefix P` writes `P.stats.json` and `P.stats.tsv`: counts by class and effect,
Ts/Tv, the substitution spectrum, indel lengths, density per contig, an allele-fraction histogram, the share
of calls with a fraction from 0.2 to 0.8, and clustering.

`qmutfind report --vcf FILE --gff GFF --region CONTIG:START-END --out page.html` writes one HTML file with
the gene models and one circle per variant, coloured by effect.

`qmutfind call --depth-out depth.bed [--depth-window 10000]` also writes the mean depth per window as BED text with four
columns (contig, start, end, depth), the same layout as a `mosdepth` regions file. It needs the built-in caller.
`qmutfind depth-plot --bed depth.bed [--vcf FILE] [--reference FA] [--order input|length|name] --out depth.svg` draws a
Manhattan-style plot of the depth along the contigs (median line, high bins clipped and marked, a tick for each call when a
VCF is given). It reads our BED or a `mosdepth` regions or per-base file, plain or `.gz`. `--order name` sorts contigs
as Chr1, Chr2, Chr10.

Each call has a `QUAL` (Phred-scaled probability that this many alternate reads come from sequencing error alone, from a
binomial model with the error rate estimated from the whole sample; it does not cover mapping errors, so a call in a
repeat can still have a high QUAL) and `DPR` (depth at the call over the sample's median depth). `--min-qual X` and
`--max-depth-ratio X` drop calls; both are off by default because they have not been tuned. With FreeBayes the QUAL is the
one FreeBayes writes.

`qmutfind plot --vcf FILE [--reference FA] [--order input|length|name] --out genome.svg` draws every call along the contigs: calls per bin on top and
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

### Loading data from the web, and where it came from

Under each input the page says where the data comes from: `Local file: name (size)`, or `From the example`, `From the
genome menu`, `From ENA run SRR...` or `From the address ...` for downloads. A file chosen in an input replaces a download.
"Load data from the web" takes a run accession (SRR, ERR or DRR; the page looks it up in ENA's lookup service, shows the
organism, platform, layout, read count and size, and downloads the FASTQ files) or https addresses of a reference, a GFF3 and
reads. Both work only when the data server allows browser downloads (NCBI and ENA do; FungiDB does not). Limits: https only,
no user name or password in an address, 2 GB per download, a file loaded whole into memory (reference, genes, VCF, coverage)
may unpack to at most 1 GB, a Cancel button, and a check of the file size and (for a run) the read count. ENA's MD5 sums
are not checked. `?limit_mb=N` and `?download_limit_mb=N` lower the two limits.

### Genome menu

`site/genomes/genomes.json` lists 155 NCBI fungal genomes (made by `bench/make_genome_list.py` from the NCBI tables; the page
downloads them from NCBI): every RefSeq Chromosome or Complete Genome assembly up to 70 Mb, 50 hand-chosen RefSeq
genomes, and two GenBank genomes, *A. fumigatus* A1163 and CEA10. CEA10 has no GFF3 at NCBI, so the page loads its
reference only. `site/genomes/fungidb.json` lists the FungiDB organisms; FungiDB blocks
downloads from web pages, so for those the page shows the `qmutfind fetch` command. `?genomes=PATH` uses another list,
but only a path on the site.
The page asks two questions. (1) Where is your reference genome? Choose the NCBI list, upload a FASTA and a GFF3, or give https
addresses for them. (2) What is your sample set? Add one row per sample with `+` and remove one with `-`. Each sample has a name
and one or more reads files from one source: an SRA run accession, https addresses (one per line), or files from your computer.
A read limit box (in thousands of reads per file, default 2,000 = 2 M; 0 uses all) keeps the first reads of each file and stops
web downloads there; the first reads are not a random sample. The samples run one after the other, each in its own worker. You
can then look at one sample (statistics, plot, gene effects, locus drawing), combine the finished samples into one table with an
optional parent sample, and export the catalog as a zip. The command line has the same combine:
`qmutfind combine --strain NAME=VCF[,BED[,EVIDENCE_VCF]] ... [--parent NAME] --out-prefix P` (`qmutfind call --evidence-out FILE`
writes the evidence VCF).

### Example data

`site/example/manifest.json` lists the examples the page can load: the Af293 genome and gene models from NCBI RefSeq
(downloaded by the browser; NCBI allows this), and reads of an evolved-strain sample (SRR6277229). The small example keeps
1 MB of reads in this repository (`site/example/reads.fastq.gz`, made by `bench/make_example_reads.sh`); the full example
downloads 137 MB from ENA (ENA allows this too). The page checks each download's size and the read count. See
`site/example/README.txt` for the sources. `?example=PATH` uses another manifest, but only a path on the site.

## Validation of the annotation

The annotation (effect class, protein change, gene and transcript) was compared with snpEff on a large, exhaustive set of variants. This shows
where the two programs agree and where they define a class differently. It is not a test against a truth set: both programs read the same GFF3,
and either one can be wrong.

**How it was done** (`bench/snpeff_validation/`; `run_validation.sh` has every command):

1. **Reference and genes.** *A. fumigatus* Af293, RefSeq `GCF_000002655.1` (FASTA and GFF3), the same files that the genome menu and `qmutfind fetch` use.
2. **snpEff database.** snpEff 5.1, built from the same FASTA and GFF3 (`build -gff3 -noCheckCds -noCheckProtein`; the CDS and protein check files do
   not exist for a GFF3-only build). snpEff 5.1 takes the transcript id of an NCBI mRNA line from `Name` (`XM_...`), but the exons name their
   parent `rna-XM_...`, so no exon finds its transcript. For snpEff only, the `rna-` and `gene-` prefixes were removed from the `Parent`
   attributes. qmutfind read the original GFF3.
3. **Variants** (`make_variants.py`, seed 17): every possible SNV at every position of 150 random protein-coding genes and 300 bp on each
   side; every possible SNV within 12 bp of each side of every exon boundary of 600 other genes (splice sites, UTR and intron edges); 4,000 random
   insertions and deletions of 1 to 9 bp in genes, and 3,000 near exon boundaries. All on the eight nuclear chromosomes. The VCF was left-aligned with
   `bcftools norm`, as the callers write indels. It has 1,498,949 records and 1,498,945 distinct variants (4 records are duplicates after normalization).
4. **Annotation.** `qmutfind annotate` (14.5 s) and `snpEff ann -canon -ud 0` (33.5 s). `-canon` makes snpEff report the canonical transcript, which
   is the longest CDS, as in qmutfind. `-ud 0` switches off the upstream and downstream labels, so a site outside genes is `intergenic_region`.
5. **Comparison** (`compare_snpeff.py`). For each variant the snpEff entry of the same transcript is used (else its first entry). Each qmutfind label has a
   list of snpEff Sequence Ontology terms that mean the same (`ACCEPT` in the script; `missense` is `missense_variant`, `nonsense` is `stop_gained`,
   `5UTR_derived` is `5_prime_UTR_variant`, and so on). A variant ends in one of three outcomes: **agree** (snpEff gives an accepted term),
   **convention** (they differ in a known definition, see below), or **differ**. The protein change is compared for variants that agree on the class: snpEff
   `p.Val29Ile` is written as `V29I`; for start-lost, stop-lost and frameshift only the reference residue and its position are compared, because the two
   programs write the rest differently. The gene and transcript are compared by name.

**Result** (the full report is `bench/snpeff_validation/results/af293_snpeff5.1.report.txt`):

| | SNVs | indels |
|---|---|---|
| variants | 1,491,951 | 6,994 |
| agree | 1,485,576 (99.573%) | 6,771 (96.812%) |
| differ by a known convention | 6,315 (0.423%) | 52 (0.743%) |
| differ, other | 60 (0.004%) | 171 (2.445%) |
| protein change equal (where both give one, and the class agrees) | 820,169 of 820,169 | 2,188 of 2,198 |
| same transcript | 1,246,611 of 1,247,238 (99.950%) | 5,793 of 5,799 (99.897%) |

Missense (587,950), synonymous (192,958), nonsense (31,211), stop-lost (4,540) and intergenic (244,713) SNVs agree in every case. The known conventions are:

- **Start codon.** qmutfind says `start_lost` for every change of the start codon. snpEff says `initiator_codon_variant` when the new codon is also a
  start codon of the table. All 1,002 cases are ATG to CTG or TTG (`M1L`). In the other 3,510 start-lost SNVs both say `start_lost`.
- **Splice region.** Near an exon edge snpEff gives only `splice_region_variant` for a position that qmutfind labels by the exons (4,798 variants: 1,533 SNVs as
  `5UTR_derived`, 3,213 SNVs as `3UTR_derived`, and 52 indels: 13, 20 and 19 `complex`). This is a different convention, not a different position.
- **Non-coding.** qmutfind says `non_coding` for features that are not protein-coding transcripts; snpEff says `intragenic_variant` or a non-coding term (567 SNVs, all in tRNA genes; the snpEff build warns "Exon's parent ... is a Gene instead of a transcript" for them).

What still differs, and why (checked against the lists of differences):

- **SNVs (60 variants): all in tRNA genes that have an intron.** snpEff gives `intragenic_variant` (48) or a UTR term (12, all in the gene `Afu5tC(GCA)3`) where qmutfind labels the
  intron (36) or the splice site (24). The snpEff build gave warnings for the tRNA genes (see above); qmutfind derives their introns from the exons.
- **Indels at exon edges (171 variants).** 71 are `splice_site` in qmutfind: snpEff gives `splice_region_variant&intron_variant` (36) or a frameshift or in-frame term together with
  `splice_region_variant` (35). 10 are `frameshift` (8) or in-frame (2) in qmutfind where snpEff puts `splice_donor_variant` first. A few others are stop-lost against an in-frame deletion (3)
  and similar. These are priority differences: both programs see the same event and rank its effects differently. 79 are `complex`, qmutfind's label for a variant that spans an exon boundary;
  snpEff gives a combination of terms, so the two cannot be compared.
- **Overlapping genes and ties (627 SNVs).** snpEff reports the same gene but another transcript, or lists it later. This happens where two genes overlap or two transcripts have the same CDS length.
  An insertion after the last base of a gene is in the gene for qmutfind and intergenic for snpEff (6 indels, all `complex` in qmutfind).
- **Protein notation.** After a frameshift snpEff 5.1 writes only the first affected residue (`p.Ala3fs`); qmutfind also gives the new residue and the stop (`A3Kfs*59`). In 10 of 2,198 frameshifts the residue or its
  position differs (by 1 to 9 residues). All 10 are near the end of the gene, and in all 10 qmutfind finds no new stop codon (`fs*?`). In 6 of them qmutfind writes `?` where the reference
  residue would be (`?471Afs*?`), which is a display defect. The cause of the different position was not examined.

**Earlier and smaller checks** (details in [PLAN.md](PLAN.md)): the effect labels were also compared with `bcftools csq` 1.25 and with snpEff 5.4c on 3,700
simulated truth variants, which found and fixed four classifier flaws (for example, the weaker label for a variant in the UTR of one gene and the CDS of another); on
real EVOL20 calls no variant disagreed. `bcftools stats` gave the same SNP, indel and Ts/Tv counts as the statistics of `qmutfind stats`.

**Limits.** One genome (Af293, standard genetic code, no mitochondrion), one snpEff version for the large run (5.1), and both programs read the same GFF3.
The GFF3 has no partial CDS with a phase other than 0, so qmutfind (which ignores the GFF phase column) was not tested on one. The variants are not real calls: they cover
every position evenly, so the percentages describe the annotation, not the calling.

**To repeat it:** `REPO=... WORK=... FASTA_GZ=... GFF_GZ=... bash bench/snpeff_validation/run_validation.sh` (needs snpEff, Java, bcftools, samtools and Python 3.12). The helper
functions of the comparison have tests in `bench/tests/test_compare_snpeff.py` (`pixi run test-bench`).

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
- `start_lost` is every change of the start codon, also to another start codon of the table (ATG to CTG or TTG); snpEff calls those `initiator_codon_variant`. For a frameshift at the stop codon the protein change shows `?` for the reference residue.
- Licence: GPL-3.0 (see [LICENSE](LICENSE)).

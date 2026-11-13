# Changelog

All notable changes to qmutfind. The project has no release yet: everything below is in development version 0.1.0. Dates are 2026-11-12 to 2026-11-13.
The measurements behind each entry are in [PLAN.md](PLAN.md).

## Unreleased

### Added: command line and core

- **`qmutfind confirm`** (core module `confirm`): confirm a knock-out or replacement allele. Reads are mapped to the genome plus the introduced cassette; split reads give the junctions where the
  cassette joins the genome, and nine checks (integration found, both sides joined, single site, at the target, target replaced, cassette coverage, cassette sequence, cassette structure, copy number) give a
  verdict (confirmed, not confirmed, inconclusive) with the evidence. Target by region or by gene name from a GFF3. Outputs: text and JSON report, junction table, depth of the cassette. Warns when the construct holds homology arms.
  See the README section "Confirming a knock-out or a replacement allele". A browser page for it is planned.
- **`qmutfind context`** and the WASM function `context`: the sequence around one variant as aligned rows (position ruler, a marker under the changed column, reference DNA,
  reference protein, variant DNA, variant protein). Exons are capitals, other bases small letters. An insertion adds a gap column, a deletion a gap in the variant row, a minus-strand
  gene is shown as it reads, a codon split by an intron is translated across the junction, and a frameshift reads on into the 3' UTR. Text or HTML (`--html`).
- **`qmutfind combine`** and WASM `combine_tables`: one table of several samples (rows are variants, one cell per sample: call, weak, ref, low, unknown, other), with an optional
  parent (status new, inherited, lost or unknown), a merged VCF of all samples (`NS`, `SAMPLES`, `AFS`, `NEW`), and a summary table.
- **Per-site evidence** for the parent status: `call --evidence-out` writes the weak signals (depth 2, fraction 0.1) as a second VCF; `--depth-out` writes windowed coverage (mosdepth-style BED) and a coverage plot.
- **`qmutfind plot`, `depth-plot`, `report`** (locus drawing with gene models, direction arrows and a variant table), and **`qmutfind fetch`** (FungiDB genomes with an API key from `EUPATHDB_KEY`).
- **Annotation:** frameshift protein change and truncation (`E3Dfs*4`, `TRUNC=new/reference`); gene names (`GENENAME`) and products (`PRODUCT`); the NCBI transcript accession from the GFF3 (`NUCCORE`, from a
  `GenBank:` entry in `Dbxref` or a `transcript_id`); the effect is judged on the canonical transcript (the longest CDS); a variant in two genes gets the stronger label; deletions are judged by the deleted bases.
- **Calling:** QUAL (binomial error model) and a depth-ratio filter (`--min-qual`, `--max-depth-ratio`, off by default); a minimum number of reads for an indel (default 6 for long reads); indels are left-aligned on the built-in path;
  `--aligner bwa-mem2` with `--caller builtin|freebayes`.
- **Pair-aware mapping** for the `rammap` path: `call --paired` and `--max-insert` (off by default; no gain was measured on real or simulated reads, see PLAN.md 12o; it stays unused and has one open defect).
- **`FastqChunker` limit:** an unfinished record may hold at most 64 MiB (an error that says the input may not be FASTQ), and a chunk is scanned once, not again at every push.

### Added: web page

- **Reference panel** with three sources: the genome menu (155 NCBI genomes with gene models where they exist, and 317 FungiDB organisms shown as a command), upload, or web addresses (a Dropbox share link is changed to the host that allows browser downloads).
- **Sample set:** rows with a name, a source (SRA run through ENA, web addresses or local files) and a note; `+` and `-` buttons; a read limit in thousands of reads per file (default 2,000); a queue that runs one sample at a time in a worker;
  the size and read-count checks of every download; an example set.
- **Combine and catalog:** combine with a parent; sortable combined table; "All samples together" (the merged VCF through summary and annotation, with a "Shared by" filter and sharing columns); "Subtract the parent"; catalog export and import (zip).
- **Saved results in the browser (IndexedDB):** finished samples, the combine and the reference are saved after each sample and each combine and restored on request. Never the reads. On by default, can be switched off, can be deleted. Results that were saved earlier and not restored are not
  overwritten (also not by a second tab); "Save now" replaces them. Without IndexedDB the page works as before.
- **Sequence under the locus drawing:** press Draw on a variant to get the aligned sequence view, with a setting for the number of bases on each side (default 40).
- **Gene links:** "Link the Gene column to" FungiDB (`https://fungidb.org/gene/GENE`) or to the NCBI nucleotide record of the transcript (only when the GFF3 gives an accession; otherwise no link).
- **Layouts:** Standard, Wide and Phone (automatic on a narrow window), five styles, and a variant table whose headers all sort.
- **Plots:** contig labels that never overlap, gene arrows and chevrons in the locus drawing, readable colours in every style (contrast-tested).
- **Build stamp** on every file that the page loads, so a cached old script cannot meet a new page.
- **Gzip files with several members** are read (fflate), and every member is now checked (CRC-32 and length).

### Added: validation and tools

- **Comparison with snpEff** on 1,498,945 distinct variants of the Af293 genome (all SNVs in 150 genes and at the exon boundaries of 600 more, and 7,000 indels): 99.57% of SNVs and 96.81% of indels agree, 0.42% and 0.74% differ by a known
  convention, 0.004% and 2.4% differ otherwise; protein changes are equal in 820,169 of 820,169 SNVs. How it was done and what differs: README, section "Validation of the annotation"; scripts in `bench/snpeff_validation/`.
- Earlier checks of the effect labels against `bcftools csq` 1.25 and snpEff 5.4c, a caller comparison on 13 samples of PRJNA1094460 (rammap, bwa-mem2 with the built-in caller, FreeBayes, bcftools), and held-out simulated truth sets (seeds 3 to 7): PLAN.md.
- A large run of the page in headless Firefox (4,000,000 reads in 150 s, 2.5 GB resident memory) and a speed comparison with the commit before the review work: PLAN.md 12s.

### Changed

- The page was redesigned (sample set, reference panel, one coloured results group). The old step 1 and the old web panel are gone.
- The WASM build uses wasm-bindgen 0.2.130.
- Reads that a read limit cuts are unpacked, counted and packed again in a worker, not on the main thread.
- A sample's files are counted together against the 2 GB download limit; with a read limit the kept reads are counted.
- The effect of a variant over two genes, an indel anchored outside a transcript and an in-frame indel that gains a stop or loses the start codon are now judged correctly (found by the `bcftools csq` check).
- Reading gzip in the browser takes about 41% longer than before because of the checksum (3.3% of the run time of a 4,000,000-read sample).

### Fixed

- A failed run freed its session only on success (about 460 MB lost per failed run); page address options are checked by origin; effect names are escaped in tables.
- Restore, import and combine are refused while a sample runs; a restored reference wins over a file chosen earlier; two tabs can no longer overwrite each other's saved results.
- A gzip file with several members failed the checksum when a chunk was larger than 64 KiB.
- Known and not fixed: for a frameshift at the stop codon the protein change shows `?` for the reference residue; pairing can promote a poor secondary alignment.

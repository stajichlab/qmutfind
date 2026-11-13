#!/usr/bin/env bash
# Call each sample against A1163 (GCA_000150145.1) with the built-in caller: VCF, coverage (1 kb windows) and evidence VCF.
# Arguments: SRR run, name. Paths are absolute (SLURM). No BASH_SOURCE is used.
set -euo pipefail
D=/bigdata/stajichlab/jstajich/projects/qmutfind_data
P=$D/prjna1094460
BIN=$D/bin/qmutfind-4734f2e
RUN=$1; NAME=$2
$BIN call --reference $P/ref/GCA_000150145.1_ASM15014v1_genomic.fna.gz \
  --reads $P/fastq/${RUN}_1.fastq.gz --reads2 $P/fastq/${RUN}_2.fastq.gz --sample $NAME --preset short \
  --out $P/calls/$NAME.vcf --depth-out $P/calls/$NAME.depth.bed --depth-window 1000 --evidence-out $P/calls/$NAME.evidence.vcf
echo "$NAME done"

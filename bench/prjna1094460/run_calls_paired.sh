#!/usr/bin/env bash
# Like run_calls.sh but with --paired. Output goes to calls_paired/. Arguments: SRR run, name. Prints the wall time and peak memory.
set -euo pipefail
D=/bigdata/stajichlab/jstajich/projects/qmutfind_data
P=$D/prjna1094460
BIN=$D/bin/qmutfind-pairing
RUN=$1; NAME=$2
mkdir -p $P/calls_paired
/usr/bin/time -f "$NAME wall_s=%e maxrss_kb=%M" $BIN call --reference $P/ref/GCA_000150145.1_ASM15014v1_genomic.fna.gz \
  --reads $P/fastq/${RUN}_1.fastq.gz --reads2 $P/fastq/${RUN}_2.fastq.gz --paired --sample $NAME --preset short \
  --out $P/calls_paired/$NAME.vcf --depth-out $P/calls_paired/$NAME.depth.bed --depth-window 1000 --evidence-out $P/calls_paired/$NAME.evidence.vcf

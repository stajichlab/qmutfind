#!/usr/bin/env bash
# Like run_calls.sh but with. Output goes to calls_single/. Arguments: SRR run, name. Prints the wall time and peak memory.
set -euo pipefail
D=/bigdata/stajichlab/jstajich/projects/qmutfind_data
P=$D/prjna1094460
BIN=$D/bin/qmutfind-pairing
RUN=$1; NAME=$2
mkdir -p $P/calls_single
/usr/bin/time -f "$NAME wall_s=%e maxrss_kb=%M" $BIN call --reference $P/ref/GCA_000150145.1_ASM15014v1_genomic.fna.gz \
  --reads $P/fastq/${RUN}_1.fastq.gz --reads2 $P/fastq/${RUN}_2.fastq.gz --sample $NAME --preset short \
  --out $P/calls_single/$NAME.vcf --depth-out $P/calls_single/$NAME.depth.bed --depth-window 1000 --evidence-out $P/calls_single/$NAME.evidence.vcf

#!/usr/bin/env bash
# First N spots of each run of PRJNA1094460 (a subsample, not a random sample), straight from NCBI with fastq-dump -X.
# Deviation from the usual prefetch route: only the first reads are needed. Each mate must have exactly N reads.
set -uo pipefail
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/sratoolkit/3.2.0/bin:$PATH
N=${N:-2000000}
OUT=/bigdata/stajichlab/jstajich/projects/qmutfind_data/prjna1094460/fastq
RUN=$1
cd "$OUT"
if [ -s ${RUN}_1.fastq.gz ] && [ -s ${RUN}_2.fastq.gz ] && [ -e ${RUN}.ok ]; then echo "$RUN already done"; exit 0; fi
rm -f ${RUN}_1.fastq.gz ${RUN}_2.fastq.gz ${RUN}.ok
fastq-dump -X $N --split-files --gzip --outdir "$OUT" $RUN > $RUN.dump.log 2>&1 || { echo "$RUN FAILED fastq-dump"; exit 1; }
c1=$(( $(zcat ${RUN}_1.fastq.gz | wc -l) / 4 )); c2=$(( $(zcat ${RUN}_2.fastq.gz | wc -l) / 4 ))
if [ "$c1" -eq "$N" ] && [ "$c2" -eq "$N" ]; then echo "$RUN ok $c1 $c2"; touch ${RUN}.ok; else echo "$RUN COUNT MISMATCH $c1 $c2 (expected $N)"; exit 1; fi

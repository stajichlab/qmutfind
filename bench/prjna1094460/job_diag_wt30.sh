#!/usr/bin/env bash
#SBATCH --job-name=diag_wt30
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/qmutfind_data/prjna1094460/diag/job_diag_wt30.%j.log
# WT30 (SRR28508724?) per-read comparison: bwa-mem2 BAM versus rammap single-read alignments. No BASH_SOURCE.
set -euo pipefail
export PATH=/bigdata/stajichlab/jstajich/projects/qmutfind/.pixi/envs/default/bin:$PATH
D=/bigdata/stajichlab/jstajich/projects/qmutfind_data
P=$D/prjna1094460
RUN=$(awk -F'\t' '$2=="WT30"{print $1}' $P/samples.tsv)
cd $P/diag
REF=$P/ref/A1163.plain.fa
bwa-mem2 index $REF 2> index.err
( bwa-mem2 mem -t 6 $REF $P/fastq/${RUN}_1.fastq.gz $P/fastq/${RUN}_2.fastq.gz 2> mem.err | samtools sort -@2 -o WT30.bwa.bam - ) &
$D/bin/map_dump $REF $P/fastq/${RUN}_1.fastq.gz 1 | gzip > WT30.rammap.m1.tsv.gz &
$D/bin/map_dump $REF $P/fastq/${RUN}_2.fastq.gz 2 | gzip > WT30.rammap.m2.tsv.gz &
wait
samtools index WT30.bwa.bam
echo done

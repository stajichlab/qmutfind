#!/usr/bin/env bash
#SBATCH --job-name=prjna1094460_compare
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/qmutfind_data/prjna1094460/compare/job_compare.%j.log
set -euo pipefail
P=/bigdata/stajichlab/jstajich/projects/qmutfind_data/prjna1094460
export PATH=/bigdata/stajichlab/jstajich/projects/qmutfind/.pixi/envs/default/bin:$PATH
zcat $P/ref/GCA_000150145.1_ASM15014v1_genomic.fna.gz > $SCRATCH/A1163.fa
bwa-mem2 index $SCRATCH/A1163.fa > /dev/null 2>&1
samtools faidx $SCRATCH/A1163.fa
# 4 samples at a time, 2 threads each: about 13 samples x 10 min / 4 = about an hour (an estimate; not yet measured)
tail -n +2 $P/samples.tsv | cut -f1,2 | xargs -P4 -L1 bash -c 'bash '$P'/compare_one.sh "$0" "$1" 2'

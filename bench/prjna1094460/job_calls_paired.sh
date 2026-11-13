#!/usr/bin/env bash
#SBATCH --job-name=prjna1094460_paired
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/qmutfind_data/prjna1094460/job_calls_paired.%j.log
set -euo pipefail
P=/bigdata/stajichlab/jstajich/projects/qmutfind_data/prjna1094460
tail -n +2 $P/samples.tsv | cut -f1,2 | xargs -P8 -L1 bash $P/run_calls_paired.sh

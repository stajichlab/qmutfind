#!/usr/bin/env bash
#SBATCH --job-name=paired_heldout
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/bigdata/stajichlab/jstajich/projects/qmutfind_data/heldout/paired_heldout.%j.log
# Simulated pairs (ART, seeds 3-5 standard truth, 10x and 20x, from the earlier comparison runs): call with and without --paired,
# same binary, and score against the truth. Rows go to heldout/paired_results.tsv with mode = single | paired. No BASH_SOURCE.
set -euo pipefail
export PATH=/bigdata/stajichlab/jstajich/projects/qmutfind/.pixi/envs/default/bin:$PATH
R=/bigdata/stajichlab/jstajich/projects/qmutfind
D=/bigdata/stajichlab/jstajich/projects/qmutfind_data
H=$D/heldout
BIN=$D/bin/qmutfind-pairing
REF_GZ=$D/ref/GCF_000002655.1_ASM265v1_genomic.fna.gz
REF_PLAIN=$D/ref/Af293.plain.fa
CONTIGS=NC_007196.1,NC_007198.1,NC_007199.1
OUT=$H/paired_results.tsv
cd $R/bench
for seed in 3 4 5; do
  for cov in 10 20; do
    C=$H/seed${seed}_cmp
    for mode in single paired; do
      flag=""; [ $mode = paired ] && flag="--paired"
      V=$H/seed${seed}_cmp/art${cov}.rammap_$mode.vcf
      $BIN call --reference $REF_GZ --reads $C/art${cov}_1.fq.gz --reads2 $C/art${cov}_2.fq.gz $flag \
        --sample s${seed}_art${cov} --preset short --out $V
      bcftools norm -f $REF_PLAIN -m -any $V -o $V.norm.vcf 2>/dev/null
      python3 score_heldout_run.py --truth $H/seed${seed}_standard/truth.norm.vcf --calls $V.norm.vcf --out $OUT \
        --seed $seed --mode $mode --tech art --coverage $cov --contigs $CONTIGS
    done
  done
done

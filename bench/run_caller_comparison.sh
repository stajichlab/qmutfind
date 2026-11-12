#!/usr/bin/env bash
# Compare callers on the held-out truth sets (standard seeds). Appends scored rows to $DATA/heldout/callers.tsv.
#   bash bench/run_caller_comparison.sh <seed>
# Run from the repository root through `pixi run`. Settings come from environment variables:
#   QMUTFIND  the qmutfind binary (default $PWD/target/release/qmutfind)
#   DATA      the data directory  (default $PWD/../qmutfind_data)
#   THREADS   threads for the external tools (default 3)
#   COVS      ART coverages (default "10 20"); ONT_COVS Nanopore coverages (default "12 24")
# Short reads: ART mates are regenerated with the same seeds as run_heldout.sh, but kept as two files (pairs).
# Callers: qmutfind built-in + bwa-mem2 (pairs), qmutfind FreeBayes + bwa-mem2 (pairs), bcftools + bwa-mem2 (pairs),
# each at the qmutfind thresholds (depth 5, SNP 0.8, indel 0.5, mapq 20) and, for bcftools, also unfiltered.
# Nanopore: minimap2 + bcftools at the long-read thresholds (SNP 0.7, indel 0.5, depth 5, indel alt reads 6).
# No BASH_SOURCE is used, so this works in SLURM jobs.
set -euo pipefail
seed="${1:?usage: run_caller_comparison.sh <seed>}"
QMUTFIND="${QMUTFIND:-$PWD/target/release/qmutfind}"
DATA="${DATA:-$PWD/../qmutfind_data}"
THREADS="${THREADS:-3}"
COVS="${COVS:-10 20}"
ONT_COVS="${ONT_COVS:-12 24}"
REF_GZ="${REF_GZ:-$DATA/ref/GCF_000002655.1_ASM265v1_genomic.fna.gz}"
REF_PLAIN="${REF_PLAIN:-${SCRATCH:?SCRATCH is not set}/map/Af293.fa}"
CONTIGS="${CONTIGS:-NC_007196.1,NC_007198.1,NC_007199.1}"
H="$DATA/heldout/seed${seed}_standard"
W="$DATA/heldout/seed${seed}_cmp"
mkdir -p "$W"
RESULTS="$DATA/heldout/callers.tsv"
TMP="$SCRATCH"
lap() { echo "[seed $seed cmp] $1 at $(date +%H:%M:%S)"; }

[ -s "$REF_PLAIN.bwt.2bit.64" ] || bwa-mem2 index "$REF_PLAIN" > /dev/null 2>&1
[ -s "$REF_PLAIN.fai" ] || samtools faidx "$REF_PLAIN"

score() {  # score <calls.vcf> <tech label> <coverage>
  bcftools norm -f "$REF_PLAIN" -m -any "$1" -o "$1.norm.vcf" 2>> "$W/norm.log"
  python3 bench/score_heldout_run.py --truth "$H/truth.norm.vcf" --calls "$1.norm.vcf" --out "$RESULTS" \
    --seed "$seed" --mode standard --tech "$2" --coverage "$3" --contigs "$CONTIGS"
}

# bcftools mpileup + call, haploid. $1 = BAM, $2 = output prefix.
bcftools_calls() {
  bcftools mpileup -f "$REF_PLAIN" -q 20 -a FORMAT/AD,FORMAT/DP --max-depth 1000 "$1" 2> "$2.mpileup.err" \
    | bcftools call --ploidy 1 -mv -Ov -o "$2.raw.vcf" 2> "$2.call.err"
}
# The thresholds of qmutfind: depth, and the allele fraction of the first alternative allele.
filter_short() {  # filter_short <in.vcf> <out.vcf>
  bcftools view -i 'FMT/DP>=5 && ((TYPE="snp" && FMT/AD[0:1]/FMT/DP>=0.8) || (TYPE!="snp" && FMT/AD[0:1]/FMT/DP>=0.5))' "$1" -o "$2"
}
filter_long() {   # the long-read rule: SNP 0.7, indel 0.5, at least 6 reads for an indel
  bcftools view -i 'FMT/DP>=5 && ((TYPE="snp" && FMT/AD[0:1]/FMT/DP>=0.7) || (TYPE!="snp" && FMT/AD[0:1]/FMT/DP>=0.5 && FMT/AD[0:1]>=6))' "$1" -o "$2"
}

for cov in $COVS; do
  art_illumina -ss HS25 -i "$H/truth.mutated.fa" -p -l 150 -f "$cov" -m 400 -s 50 -rs $((seed * 100 + cov)) -na -o "$W/art${cov}_" > /dev/null 2>&1
  gzip -f "$W/art${cov}_1.fq" "$W/art${cov}_2.fq"
  r1="$W/art${cov}_1.fq.gz"; r2="$W/art${cov}_2.fq.gz"
  lap "ART ${cov}x reads ready"
  for caller in builtin freebayes; do
    "$QMUTFIND" call --reference "$REF_GZ" --reads "$r1" --reads2 "$r2" --sample "s${seed}" --preset short \
      --aligner bwa-mem2 --caller "$caller" --threads "$THREADS" --tmpdir "$TMP" --out "$W/art${cov}.bwa_$caller.vcf" 2> "$W/art${cov}.bwa_$caller.err"
    score "$W/art${cov}.bwa_$caller.vcf" "art_bwa_$caller" "$cov"
    lap "ART ${cov}x bwa-mem2 + $caller done"
  done
  bwa-mem2 mem -t "$THREADS" -R "@RG\tID:s\tSM:s" "$REF_PLAIN" "$r1" "$r2" 2> "$W/art${cov}.bwa.err" \
    | samtools sort -@ 2 -T "$TMP/sort_$seed" -o "$W/art${cov}.bam" - 2> "$W/art${cov}.sort.err"
  samtools index "$W/art${cov}.bam"
  bcftools_calls "$W/art${cov}.bam" "$W/art${cov}.bcftools"
  score "$W/art${cov}.bcftools.raw.vcf" "art_bwa_bcftools_raw" "$cov"
  filter_short "$W/art${cov}.bcftools.raw.vcf" "$W/art${cov}.bcftools.filtered.vcf"
  score "$W/art${cov}.bcftools.filtered.vcf" "art_bwa_bcftools_filtered" "$cov"
  rm -f "$W/art${cov}.bam" "$W/art${cov}.bam.bai"
  lap "ART ${cov}x bcftools done"
done

for cov in $ONT_COVS; do
  minimap2 -ax map-ont -t "$THREADS" "$REF_PLAIN" "$H/ont${cov}.fastq.gz" 2> "$W/ont${cov}.mm2.err" \
    | samtools sort -@ 2 -T "$TMP/sort_$seed" -o "$W/ont${cov}.bam" - 2> "$W/ont${cov}.sort.err"
  samtools index "$W/ont${cov}.bam"
  bcftools_calls "$W/ont${cov}.bam" "$W/ont${cov}.bcftools"
  score "$W/ont${cov}.bcftools.raw.vcf" "ont_mm2_bcftools_raw" "$cov"
  filter_long "$W/ont${cov}.bcftools.raw.vcf" "$W/ont${cov}.bcftools.filtered.vcf"
  score "$W/ont${cov}.bcftools.filtered.vcf" "ont_mm2_bcftools_filtered" "$cov"
  rm -f "$W/ont${cov}.bam" "$W/ont${cov}.bam.bai"
  lap "ONT ${cov}x minimap2 + bcftools done"
done
lap "all done"

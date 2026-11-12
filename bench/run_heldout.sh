#!/usr/bin/env bash
# Run one held-out seed from start to finish and append the scored rows to results.tsv.
#   bash bench/run_heldout.sh <seed> <standard|hard>
# Run from the repository root through `pixi run`. Settings come from environment variables:
#   QMUTFIND   the release binary          (default $PWD/target/release/qmutfind)
#   DATA       the data directory          (default $PWD/../qmutfind_data)
#   REF_GZ     the reference, gzip FASTA   (default $DATA/ref/GCF_000002655.1_ASM265v1_genomic.fna.gz)
#   REF_PLAIN  the reference, plain FASTA  (default $SCRATCH/map/Af293.fa; used by bcftools norm)
#   CONTIGS    the contigs the reads come from, comma-separated (default chromosomes 3, 5 and 6)
# No BASH_SOURCE is used, so this works in SLURM jobs.
set -euo pipefail
seed="${1:?usage: run_heldout.sh <seed> <standard|hard>}"
mode="${2:?usage: run_heldout.sh <seed> <standard|hard>}"
QMUTFIND="${QMUTFIND:-$PWD/target/release/qmutfind}"
DATA="${DATA:-$PWD/../qmutfind_data}"
REF_GZ="${REF_GZ:-$DATA/ref/GCF_000002655.1_ASM265v1_genomic.fna.gz}"
REF_PLAIN="${REF_PLAIN:-${SCRATCH:?SCRATCH is not set}/map/Af293.fa}"
CONTIGS="${CONTIGS:-NC_007196.1,NC_007198.1,NC_007199.1}"
OUT="$DATA/heldout"
W="$OUT/seed${seed}_${mode}"
mkdir -p "$W"
RESULTS="$OUT/results.tsv"
lap() { echo "[seed $seed $mode] $1 at $(date +%H:%M:%S)"; }

# 1. The subset reference: only the contigs the reads come from.
python3 - "$REF_PLAIN" "$CONTIGS" "$W/sub.fa" <<'PY'
import sys
keep = set(sys.argv[2].split(","))
out, on = open(sys.argv[3], "w"), False
for line in open(sys.argv[1]):
    if line.startswith(">"):
        on = line[1:].split()[0] in keep
    if on:
        out.write(line)
PY

# 2. The truth set on the subset, then normalized against the whole reference.
if [ "$mode" = hard ]; then
  hard=(--min-spacing 30 --homopolymer-fraction 0.5 --cluster-fraction 0.25)
else
  hard=()
fi
python3 bench/make_truth.py --reference "$W/sub.fa" --out-prefix "$W/truth" --seed "$seed" \
  --snps 1500 --small-indels 400 --large-indels 100 "${hard[@]}"
bcftools norm -f "$REF_PLAIN" -m -any "$W/truth.truth.vcf" -o "$W/truth.norm.vcf" 2> "$W/truth.norm.log"
lap "truth done"

score() {  # score <calls.vcf> <tech> <coverage>
  bcftools norm -f "$REF_PLAIN" -m -any "$1" -o "$1.norm.vcf" 2>> "$W/norm.log"
  python3 bench/score_heldout_run.py --truth "$W/truth.norm.vcf" --calls "$1.norm.vcf" --out "$RESULTS" \
    --seed "$seed" --mode "$mode" --tech "$2" --coverage "$3" --contigs "$CONTIGS"
}

# 3. Short reads: ART at four coverages, each its own simulation, both mates mapped as single-end.
for cov in 5 10 20 30; do
  art_illumina -ss HS25 -i "$W/truth.mutated.fa" -p -l 150 -f "$cov" -m 400 -s 50 -rs $((seed * 100 + cov)) -na -o "$W/art${cov}_" > /dev/null 2>&1
  cat "$W/art${cov}_1.fq" "$W/art${cov}_2.fq" | gzip -c > "$W/art${cov}.fq.gz"
  rm -f "$W/art${cov}_1.fq" "$W/art${cov}_2.fq"
  "$QMUTFIND" call --reference "$REF_GZ" --reads "$W/art${cov}.fq.gz" --sample "s${seed}_art${cov}" --preset short --out "$W/art${cov}.vcf" 2> "$W/art${cov}.err"
  score "$W/art${cov}.vcf" art "$cov"
  lap "ART ${cov}x done"
done

# 4. Nanopore reads: four independent 6x simulations. 6x is one of them, 12x two, 24x all four.
for i in 1 2 3 4; do
  badread simulate --reference "$W/truth.mutated.fa" --quantity 6x --length 8000,6000 --seed $((seed * 10 + i)) > "$W/ont_$i.fastq" 2> "$W/ont_$i.err" &
done
wait
lap "Badread done"
declare -A set_for=( [6]="1" [12]="1 2" [24]="1 2 3 4" )
for cov in 6 12 24; do
  : > "$W/ont${cov}.fastq"
  for i in ${set_for[$cov]}; do cat "$W/ont_$i.fastq" >> "$W/ont${cov}.fastq"; done
  gzip -f "$W/ont${cov}.fastq"
  "$QMUTFIND" call --reference "$REF_GZ" --reads "$W/ont${cov}.fastq.gz" --sample "s${seed}_ont${cov}" --preset long --out "$W/ont${cov}.vcf" 2> "$W/ont${cov}.err"
  score "$W/ont${cov}.vcf" ont "$cov"
  lap "ONT ${cov}x done"
done
rm -f "$W"/ont_?.fastq
lap "all done"

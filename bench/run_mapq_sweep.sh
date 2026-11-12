#!/usr/bin/env bash
# Sweep --min-mapq for the bwa-mem2 paths (built-in caller and FreeBayes). Writes one VCF per run.
#   bash bench/run_mapq_sweep.sh
# Run from the repository root through `pixi run`. Settings come from environment variables:
#   QMUTFIND  the qmutfind binary   (default $PWD/target/release/qmutfind)
#   DATA      the data directory    (default $PWD/../qmutfind_data)
#   THREADS   threads for bwa-mem2  (default 2)
#   QS        mapq values            (default "0 10 20 30 40")
# Datasets: the real parent (SRR6277228) and EVOL20_3 (SRR6277229), and the simulated ART seed-2 pairs.
# Note: for FreeBayes, q0 is sent as 1 (`--min-mapping-quality` 1), so q0 and q1 are the same run.
# No BASH_SOURCE is used, so this works in SLURM jobs.
set -euo pipefail
QMUTFIND="${QMUTFIND:-$PWD/target/release/qmutfind}"
DATA="${DATA:-$PWD/../qmutfind_data}"
THREADS="${THREADS:-2}"
QS="${QS:-0 10 20 30 40}"
REF="$DATA/ref/GCF_000002655.1_ASM265v1_genomic.fna.gz"
OUT="$DATA/mapq_bwa"
mkdir -p "$OUT"
TMP="${SCRATCH:?SCRATCH is not set}"

run() {  # run <name> <reads1> <reads2> <sample> <caller> <q>
  local name="$1" r1="$2" r2="$3" sample="$4" caller="$5" q="$6"
  local out="$OUT/$name.bwa_$caller.q$q.vcf"
  if [ -s "$out" ]; then echo "skip $out"; return; fi
  local start=$SECONDS
  "$QMUTFIND" call --reference "$REF" --reads "$r1" --reads2 "$r2" --sample "$sample" --preset short \
    --aligner bwa-mem2 --caller "$caller" --min-mapq "$q" --threads "$THREADS" --tmpdir "$TMP" \
    --out "$out.tmp" 2> "$OUT/$name.bwa_$caller.q$q.err"
  mv "$out.tmp" "$out"
  echo "$name $caller q$q: $(grep -vc '^#' "$out") calls, $((SECONDS - start)) s"
}

for q in $QS; do
  for caller in builtin freebayes; do
    run SRR6277228 "$DATA/fastq/SRR6277228_1.fastq.gz" "$DATA/fastq/SRR6277228_2.fastq.gz" parent "$caller" "$q"
    run SRR6277229 "$DATA/fastq/SRR6277229_1.fastq.gz" "$DATA/fastq/SRR6277229_2.fastq.gz" EVOL20_3 "$caller" "$q"
    run af293_s2_art "$DATA/bench/af293_s2_art_1.fq.gz" "$DATA/bench/af293_s2_art_2.fq.gz" art "$caller" "$q"
  done
done

#!/usr/bin/env bash
# Four more call sets for one PRJNA1094460 sample, against A1163, with bwa-mem2 mapping the two mates as a pair:
#   builtin  = qmutfind built-in caller on the bwa-mem2 alignments      (CLI --aligner bwa-mem2 --caller builtin)
#   freebayes = qmutfind FreeBayes path on the same mapping             (CLI --caller freebayes)
#   bcftools raw and filtered: bwa-mem2 | samtools sort | bcftools mpileup | bcftools call --ploidy 1 -mv, then the
#   qmutfind thresholds (depth 5, SNP fraction 0.8, indel fraction 0.5).
# Arguments: SRR run, sample name, threads. Absolute paths; no BASH_SOURCE (SLURM).
set -euo pipefail
RUN=$1; NAME=$2; T=${3:-2}
D=/bigdata/stajichlab/jstajich/projects/qmutfind_data
P=$D/prjna1094460; O=$P/compare
BIN=$D/bin/qmutfind-4734f2e
export PATH=/bigdata/stajichlab/jstajich/projects/qmutfind/.pixi/envs/default/bin:$PATH
REF_GZ=$P/ref/GCA_000150145.1_ASM15014v1_genomic.fna.gz
W=${SCRATCH:?SCRATCH is not set}/cmp_$NAME; mkdir -p $W
REF=$SCRATCH/A1163.fa
r1=$P/fastq/${RUN}_1.fastq.gz; r2=$P/fastq/${RUN}_2.fastq.gz
lap() { echo "[$NAME] $1 at $(date +%H:%M:%S)"; }

for caller in builtin freebayes; do
  $BIN call --reference $REF_GZ --reads $r1 --reads2 $r2 --sample $NAME --preset short --aligner bwa-mem2 --caller $caller \
    --threads $T --tmpdir $W --out $O/$NAME.bwa_$caller.vcf 2> $O/$NAME.bwa_$caller.err
  lap "bwa-mem2 + $caller done"
done
bwa-mem2 mem -t $T -R "@RG\tID:$NAME\tSM:$NAME" $REF $r1 $r2 2> $O/$NAME.bwa.err | samtools sort -@ 2 -T $W/sort -o $W/$NAME.bam - 2> $O/$NAME.sort.err
samtools index $W/$NAME.bam
bcftools mpileup -f $REF -q 20 -a FORMAT/AD,FORMAT/DP --max-depth 1000 $W/$NAME.bam 2> $O/$NAME.mpileup.err | bcftools call --ploidy 1 -mv -Ov -o $O/$NAME.bcftools.raw.vcf 2> $O/$NAME.call.err
bcftools view -i 'FMT/DP>=5 && ((TYPE="snp" && FMT/AD[0:1]/FMT/DP>=0.8) || (TYPE!="snp" && FMT/AD[0:1]/FMT/DP>=0.5))' $O/$NAME.bcftools.raw.vcf -o $O/$NAME.bcftools.filtered.vcf
lap "bcftools done"
rm -rf $W

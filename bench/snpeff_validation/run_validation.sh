#!/usr/bin/env bash
# Compare `qmutfind annotate` with snpEff on about 1.5 million variants of the Af293 RefSeq genome.
# Usage: REPO=/path/to/qmutfind WORK=/path/to/workdir bash run_validation.sh
# Needs: the qmutfind CLI (REPO/target/release/qmutfind, or set QMUTFIND), snpEff 5.1 (on UCR HPCC: `module load snpEff/5.1`, which sets
# SNPEFFJAR and SNPEFFDIR; elsewhere set both), Java, bcftools, samtools and Python 3.12. The reference and GFF3 are
# GCF_000002655.1 (the files that `qmutfind fetch` or the genome menu use). No BASH_SOURCE is used, so it also runs in SLURM jobs.
set -euo pipefail
REPO="${REPO:?set REPO to the qmutfind repository}"
WORK="${WORK:?set WORK to an empty working directory}"
QMUTFIND="${QMUTFIND:-$REPO/target/release/qmutfind}"
FASTA_GZ="${FASTA_GZ:?set FASTA_GZ to GCF_000002655.1_ASM265v1_genomic.fna.gz}"
GFF_GZ="${GFF_GZ:?set GFF_GZ to GCF_000002655.1_ASM265v1_genomic.gff.gz}"
: "${SNPEFFJAR:?load snpEff (module load snpEff/5.1) or set SNPEFFJAR and SNPEFFDIR}"
: "${SNPEFFDIR:?set SNPEFFDIR to the folder that holds snpEff.config}"
SEED="${SEED:-17}"
mkdir -p "$WORK/data/af293" && cd "$WORK"

# 1. The reference and the GFF3 as plain files. snpEff takes the transcript id of an NCBI mRNA line from its Name (XM_...), while the exons name
#    their parent as Parent=rna-XM_..., so no exon finds its transcript. Removing the `rna-` and `gene-` prefixes from the Parent
#    attributes makes the two agree. The copy for qmutfind keeps the original GFF3.
zcat "$FASTA_GZ" > data/af293/sequences.fa
zcat "$GFF_GZ" > genes.ncbi.gff
sed -E 's/(^|;)(ID|Parent)=(rna|gene)-/\1\2=/g' genes.ncbi.gff > data/af293/genes.gff
samtools faidx data/af293/sequences.fa

# 2. The snpEff database. The CDS and protein check files are not available for a GFF3-only build, so the checks are off (-noCheckCds -noCheckProtein).
cp "$SNPEFFDIR/snpEff.config" snpEff.config
sed -i "s|^data.dir.*|data.dir = $WORK/data/|" snpEff.config
echo "af293.genome : Aspergillus_fumigatus_Af293" >> snpEff.config
java -Xmx8g -jar "$SNPEFFJAR" build -c snpEff.config -gff3 -noCheckCds -noCheckProtein -v af293 > build.log 2>&1

# 3. The variants: every SNV in 150 random genes and 300 bp around them, every SNV near the exon boundaries of 600 more genes, and 7,000 random indels.
#    Left-aligned with bcftools norm, as the callers write them.
python3.12 "$REPO/bench/snpeff_validation/make_variants.py" genes.ncbi.gff data/af293/sequences.fa "$SEED" raw.vcf
bcftools norm -f data/af293/sequences.fa -m -any raw.vcf -o variants.vcf 2> norm.log

# 4. Both annotators. -canon: snpEff's canonical transcript is the longest CDS, as in qmutfind. -ud 0: no upstream or downstream labels.
"$QMUTFIND" annotate --vcf variants.vcf --reference "$FASTA_GZ" --gff "$GFF_GZ" --out qmutfind.vcf
java -Xmx12g -jar "$SNPEFFJAR" ann -c snpEff.config -noStats -canon -ud 0 af293 variants.vcf > snpeff.vcf 2> snpeff.err

# 5. The comparison: report and the lists of differences.
python3.12 "$REPO/bench/snpeff_validation/compare_snpeff.py" qmutfind.vcf snpeff.vcf cmp
echo "report: $WORK/cmp.report.txt"

#!/usr/bin/env bash
# Make the small example read set of the site: the reads of a sample that map within WINDOW bases of each locus in a
# list, taken from a BAM, as one gzip FASTQ file.
#   bash bench/make_example_reads.sh <sample.bam> <loci.tsv> <out.fastq.gz> [window]
# loci.tsv has two columns, contig and 1-based position. Windows that touch are merged, so a read is written once.
# Both mates of a pair are written when they map in a window (the page maps reads one by one, as single-end reads).
# Run through `pixi run` (needs samtools). No BASH_SOURCE is used.
set -euo pipefail
bam="${1:?usage: make_example_reads.sh <bam> <loci.tsv> <out.fastq.gz> [window]}"
loci="${2:?}"
out="${3:?}"
window="${4:-25000}"
regions="${SCRATCH:?SCRATCH is not set}/example_regions.txt"
sort -k1,1 -k2,2n "$loci" | awk -v w="$window" '
  { s = $2 - w; if (s < 1) s = 1; e = $2 + w
    if ($1 == c && s <= pe + 1) { if (e > pe) pe = e } else { if (c != "") print c ":" ps "-" pe; c = $1; ps = s; pe = e } }
  END { if (c != "") print c ":" ps "-" pe }' > "$regions"
# -F 0x904: skip unmapped, secondary and supplementary records.
samtools view -b -F 0x904 "$bam" $(cat "$regions") | samtools fastq -n - 2> /dev/null | gzip -9 > "$out"
echo "regions: $(wc -l < "$regions"); reads: $(( $(gzip -dc "$out" | wc -l) / 4 ))" >&2

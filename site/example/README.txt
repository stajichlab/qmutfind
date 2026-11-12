Example data of the qmutfind site
=================================

manifest.json lists the files. The page downloads them in your browser when you press a button.

Reference and genes (both examples)
  NCBI RefSeq GCF_000002655.1 (Aspergillus fumigatus Af293), genomic FASTA and GFF3, from
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/655/GCF_000002655.1_ASM265v1/
  NCBI allows downloads from web pages (Access-Control-Allow-Origin: *), so nothing is copied into this repository.

Full example reads
  ENA run SRR6277229 (BioProject PRJNA417720, sample EVOL20_3), paired Illumina reads, 1,427,487 read pairs.
  https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR627/009/SRR6277229/SRR6277229_{1,2}.fastq.gz
  ENA allows downloads from web pages. The two files are 63,389,773 and 73,227,390 bytes. The page checks the sizes
  and, after mapping, the read count (2,854,974 reads; the page maps each mate as a separate read).

Small example reads (reads.fastq.gz, kept in this repository)
  The reads of SRR6277229 whose alignment lies within 25 kb of 8 loci, as one gzip FASTQ file (33,729 reads). The loci
  are the variants that the first run of qmutfind found in this sample and not in its parent sample SRR6277228.
  Made with: bash bench/make_example_reads.sh SRR6277229.bam loci.tsv reads.fastq.gz 25000
  where the BAM is a minimap2 short-read alignment of both mates to the RefSeq genome.
  Only these regions have reads. Everywhere else the depth is 0, and the coverage plot shows peaks at the 8 regions.

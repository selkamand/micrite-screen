#!/usr/bin/env bash

#  We use --exclude-flags argument to get rid of alignment records that are any of the following
#   1. supplementary (split reads)
#   2. PCR or optical duplicates
#   3. Read fails platform/vendor quality checks
#   4. not primary alignment (secondary, i.e. multimapped)

# We use the --include-flags argument to include only reads that have at least one of the following flags
#   1. Read unmapped
#   2. Mate unmapped
threads=0

# We write our fastqs to unmapped.R1.f1, unmapped.R2.fq and unmapped.singleton.fq
samtools fastq \
  --exclude-flags 3840 \
  --include-flags 12 \
  -1 unmapped.R1.fq \
  -2 unmapped.R2.fq \
  -s unmapped.singleton.fq \
  -@ ${threads} \
  ${bam}


# Also extract reads whose reference seq is a decoy contig EXCLUDING 
# those with any of the following properties:
# 1. Read unmapped
# 2. Not primary alignment (secondary, i.e. multimapped)
# 3. Read is PCR or optical duplicate
# 4. Supplementary (split reads)

# Note we add a --fetch-pairs statement to ensure that a read pair where one read maps to
# decoy but the other maps to decoy contig (imaging reads from fragments spanning a viral integration site)
# are kept. 
samtools view --fetch-pairs \
  --exclude-flags 3844 \
  ${bam} \
  ${microbial_contigs} | \
  samtools fastq \
  -1 mapped_to_decoy.R1.fq \
  -2 mapped_to_decoy.R2.fq \
  -s mapped_to_decoy.singleton.fq \
  -@ ${threads}

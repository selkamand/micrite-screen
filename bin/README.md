# Fetching Unmapped reads


The fetch_unmapped_reads.sh script extracts unmapped reads for downstream metagenomics classification.


It generates a **single set of FASTQ files** from an aligned BAM file containing all read pairs or fragments that are likely to represent non-host or anomalous sequence and are therefore appropriate for downstream metagenomic classification.

The output includes:
- read pairs where **either read is unmapped**, and
- read pairs where **either read maps to a user-specified decoy contig** (for example microbial, viral, or synthetic references).

Only primary, high-quality alignments are considered; secondary alignments (multimaps), duplicates, QC-failed reads, and supplementary (split/chimeric) alignments are excluded. If only one read in pair fails QC you may end up with singletons.

### Considerations for integrated viruses

If one read in a pair maps to a decoy contig, the **entire pair is retained**, even if the mate maps to host or is unmapped. This is to ensure that read pairs spanning microbial contigs (where one represents mostly host and the other integrated virus) are included in the output. While supplementary alignments (representing split reads) are always excluded - we expect most split-reads will still be included since the mate pair will usually map to different species when decoys are included or include one unmapped read if no decoy is included.

### Singletons

one mate is marked QC-fail / duplicate / secondary / supplementary and gets excluded

### Assumptions

As per the SAM Format specification: 
> reads/segments having identical QNAME are regarded to come from the same template.
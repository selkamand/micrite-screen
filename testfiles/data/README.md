# README

## Reference Genome
Reference `ref/ref.fa` was manually created.

Reference Indexes were built with bwa-mem2 

```
bwa-mem2 index ref/ref.fa
```

## Read simulation and alignment

Reads paired1.fq & paired2.fq were simulated from ref.fa using ART:
```
art_illumina \
    -ss HS25 \          ## Sequencing system: HiSeq 2500
    -i ref/ref.fa \     ## Input reference FASTA
    -p \                ## Generate paired-end reads
    -l 50 \             ## Read length (bases)
    -f 20 \             ## Fold coverage (per fragment)
    -m 200 \            ## Mean fragment size
    -s 10 \             ## Std deviation of fragment size
    --noALN \           ## Do not produce an alignment file
    -o paired           ## Output prefix (paired1.fq / paired2.fq)
    

art_illumina -ss HS25 -i ref/ref.fa --paired --len 50 --fcov 20 --mflen 200 --noALN --sdev 10 -o paired
```

Simulated reads were aligned back to ref with bwa-mem2

```
bwa-mem2 mem ref/ref.fa paired1.fq paired2.fq | samtools sort -O bam > aligned.bam
samtools index aligned.bam
```

Computed stats by running `samtools idxstats aligned.bam`

This gives us a bam with

- 40 reads mapped to *refseq1*
- 20 reads mapped to *refseq2* 
- 20 reads mapped to *microbedecoy*


We then drop the decoy sequence from the reference but keep the 20 reads simulated from that sequence.
After realignment we get

```
bwa-mem2 mem ref_no_decoy/ref_no_decoy.fa paired1.fq paired2.fq | samtools sort -O bam > aligned.no_decoy.bam
samtools index aligned.no_decoy.bam
```

`aligned.no_decoy.bam` has:

- 20 unmapped reads (*)
- 60 aligned reads


### one of pair unmapped

pairsplit fastqs have 1 read pair. The first read maps to refseq1 and the second maps to refseq2
paired1.pairsplit.fq  paired2.pairsplit.fq
```
bwa-mem2 mem ref_no_decoy/ref_no_decoy.fa paired1.pairsplit.fq paired2.pairsplit.fq | samtools sort -O bam > aligned.pairsplit.bam
samtools index aligned.pairsplit.bam
```

## Databases

We supply a krakenuniq database (`minuku`) that includes all 3 sequences in ref.

We use a short minimizer-len to ensure datbase stays tiny

```
krakenuniq-build --db miniku --build --minimizer-len 10
```


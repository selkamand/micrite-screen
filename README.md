# micrite-screen

A nextflow workflow for detection of microbial reads from cancer whole-genome-sequenced bams.

Based on the Ge et al. 2025 manuscript, modified to work from oncoanalyser bam files.


## Pipeline

We start from an oncoanalyser aligned bam, extract all unmapped reads, realign to CHM13 with bowtie, pull unmapped reads once more and run kraken-uniq to classify reads against the microbial2023 database.


## Running the pipeline

### Step 1: Download required reference genomes/databases

Start by downloading the reference genomes and databases this pipeline will use.

#### i) CHM13v2.0 reference genome

Download and prepare bowtie indices for CHM13v2.0 genome.
```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz

gzip -dc chm13v2.0.fa.gz | bgzip -c > chm13v2.0.fa.bgz

bowtie2-build chm13v2.0.fa.bgz
```
#### ii) Microbial2023 database

Download microbial2023 krakenuniq database

```
wget "https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/kuniq_microbialdb_minus_kdb.20230808.tgz"
wget "https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/database.kdb"
```

Untar the database & move `database.kdb` file into database folder

```
tar -xvf kuniq_microbialdb_minus_kdb.20230808.tgz

mv database.kdb kuniq_microbialdb_minus_kdb.20230808
```

### Step 2: Run nextflow pipeline

```

nextflow run selkamand/micrite-screen \
  --bam sample.bam \
  --decoys contigs.txt \
  --ref /path/to/bowtie2/index/prefix/chm13v2.0.fa \
  --kraken_db /path/to/krakenuniq_db \
  --threads 8 \
  --threads_kraken 4 \
  --bowtie2_preset sensitive \
  --preload_size 4G \
  --outdir results
```


## Output Files

All files are written to `--outdir`.

**Pre-alignment unmapped reads**
- `<sample>.pre_R1.fastq.gz`
- `<sample>.pre_R2.fastq.gz`
  Unmapped read pairs extracted from the original BAM.

**CHM13 re-alignment**
- `<sample>.sorted.bam`
- `<sample>.sorted.bam.bai`
  Coordinate-sorted BAM and index after re-alignment to CHM13.

**Post-alignment unmapped reads**
- `<sample>.post_R1.fastq.gz`
- `<sample>.post_R2.fastq.gz`
  Read pairs still unmapped after CHM13 alignment (input to KrakenUniq).

**Taxonomic classification**
- `<sample>.krakenuniq.report.txt`
  KrakenUniq taxonomic classification report for final unmapped reads.

# Building the dockerfile

From inside this directory.

Build local version for OSX

```{bash}
docker buildx build --platform linux/arm64 --load --tag selkamandcci/micrite-screen:0.0.1 .
```

Build final version to push to dockerhub

```{bash}
docker buildx build --push --platform linux/amd64,linux/arm64 --tag selkamandcci/micrite-screen:0.0.1 .
```

# Testing 

This repo includes some small files for testing. If you git clone the repo you should be able to run the following.

You may need to change `--profile` from singularity to docker if running on MacOS

```
nextflow run . --profile singularity \
  --ref testfiles/data/ref/ref.fa \
  --kraken_db testfiles/data/dbs/minuku/ \
  --bam testfiles/data/aligned.no_decoy.bam \ 
  --decoys testfiles/data/decoys.txt \
  --outdir outdir 
```

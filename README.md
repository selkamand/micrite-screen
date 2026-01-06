# goldmicrobe

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

TODO: add after docker build
```
nextflow run -profile singularity selkamand/goldmicrobe \
 --bam=/path/to/bam \
 --database=/path/to/microbial2023 \
 --genome=/path/to/chm13_genome \
 --cores=1 \ # nthreads
 --outdir=/path/to/write/results
```


# Building the dockerfile

From inside this directory.

Build local version for OSX

```{bash}
docker buildx build --platform linux/arm64 --load --tag selkamandcci/goldmicrobe:0.0.1 .
```

Build final version to push to dockerhub

```{bash}
docker buildx build --push --platform linux/amd64,linux/arm64 --tag selkamandcci/goldmicrobe:0.0.1 .
```
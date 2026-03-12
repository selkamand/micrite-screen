nextflow.enable.dsl = 2

process FETCH_UNMAPPED_PRE {
    input:
    tuple val(sampleid), path(bam), path(bai), path(decoys)

    output:
    tuple val(sampleid), path("${sampleid}.unmapped.R1.fq.gz"), path("${sampleid}.unmapped.R2.fq.gz")

    script:
    """
    ls;
     fetch_unmapped_reads.sh \
    --bam ${bam} \
    --decoys ${decoys} \
    --prefix "${sampleid}.unmapped" \
    --threads ${task.cpus} \
    --outdir .
  """
}

process FETCH_UNMAPPED_POST {

    tag "Unmap from ${bam}"

    input:
    tuple val(sampleid), path(bam), path(bai), path(decoys)

    output:
    tuple val(sampleid), path("${sampleid}.hostdepleted.R1.fq.gz"), path("${sampleid}.hostdepleted.R2.fq.gz")

    script:
    """
    fetch_unmapped_reads.sh \
    --bam ${bam} \
    --decoys ${decoys} \
    --prefix ${sampleid}.hostdepleted \
    --threads ${task.cpus} \
    --outdir .
  """
}

process ALIGN_BOWTIE2 {
    input:
    tuple val(sampleid), path(ref), path(r1), path(r2), path(bowtie_indexes), val(preset)

    output:
    tuple val(sampleid), path("${sampleid}.realigned.sorted.bam"), path("${sampleid}.realigned.sorted.bam.bai")

    script:
    """
  align_bowtie2.sh \
    -x ${ref} \
    -1 ${r1} \
    -2 ${r2} \
    -o ${sampleid}.realigned \
    -t ${task.cpus} \
    --preset ${preset}
  """
}

process HOST_DEPLETION_STATS {
    input:
    tuple val(sampleid), path(bam_original), path(bai_original), path(decoys), path(r1_unmapped), path(r2_unmapped), path(r1_depleted), path(r2_depleted)

    output:
    tuple val(sampleid), path("${sampleid}.original.bam.stats.tsv"), path("${sampleid}.unmapped.stats.tsv"), path("${sampleid}.host_depleted.stats.tsv")

    script:
    """
    set -euo pipefail

    samtools flagstat -@ ${task.cpus} ${bam_original} > ${sampleid}.original.bam.stats.tsv

    seqkit stats -j ${task.cpus} ${r1_unmapped} ${r2_unmapped} > ${sampleid}.unmapped.stats.tsv

    seqkit stats -j ${task.cpus} ${r1_depleted} ${r2_depleted} > ${sampleid}.host_depleted.stats.tsv
    """
}

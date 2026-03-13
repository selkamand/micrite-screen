nextflow.enable.dsl = 2

process FETCH_UNMAPPED {

    tag "Unmap from ${bam}"

    input:
    val label
    tuple val(sampleid), path(bam), path(bai)
    path decoys

    output:
    tuple val(sampleid), path("${sampleid}.${label}.R1.fq.gz"), path("${sampleid}.${label}.R2.fq.gz")

    script:
    """
    fetch_unmapped_reads.sh \
    --bam ${bam} \
    --decoys ${decoys} \
    --prefix ${sampleid}.${label} \
    --threads ${task.cpus} \
    --outdir .
  """
}
process ALIGN_BOWTIE2 {
    input:
    val label
    tuple val(sampleid), path(r1), path(r2)
    tuple path(reference_fasta), path(bowtie_indexes), val(preset)

    output:
    tuple val(sampleid), path("${sampleid}.${label}.sorted.bam"), path("${sampleid}.${label}.sorted.bam.bai")

    script:
    """
  align_bowtie2.sh \
    -x ${reference_fasta} \
    -1 ${r1} \
    -2 ${r2} \
    -o ${sampleid}.${label} \
    -t ${task.cpus} \
    --preset ${preset}
  """
}


process HOST_DEPLETION_STATS_BAM {
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

process HOST_DEPLETION_STATS_FASTQ {
    input:
    tuple val(sampleid), path(r1_original), path(r2_original), path(r1_depleted), path(r2_depleted)

    output:
    tuple val(sampleid), path("${sampleid}.original.stats.tsv"), path("${sampleid}.host_depleted.stats.tsv")

    script:
    """
    set -euo pipefail

    seqkit stats -j ${task.cpus} ${r1_original} ${r2_original} > ${sampleid}.original.stats.tsv
    seqkit stats -j ${task.cpus} ${r1_depleted} ${r2_depleted} > ${sampleid}.host_depleted.stats.tsv
    """
}

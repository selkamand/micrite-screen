nextflow.enable.dsl = 2

process BAMSTATS {
    input:
    val label
    tuple val(sampleid), path(bam), path(bai)

    output:
    tuple val(sampleid), path("${sampleid}.${label}.bam.flagstats.tsv"), path("${sampleid}.${label}.bam.idxstats.tsv")

    script:
    """
    set -euo pipefail

    samtools flagstat -@ ${task.cpus} ${bam} > ${sampleid}.${label}.bam.flagstats.tsv
    samtools idxstats -@ ${task.cpus} ${bam} > ${sampleid}.${label}.bam.idxstats.tsv
    """
}

process SEQKIT_STATS {
    input:
    val label
    tuple val(sampleid), path(r1), path(r2)

    output:
    tuple val(sampleid), path("${sampleid}.${label}.stats.tsv"), path("${sampleid}.${label}.stats.tsv")

    script:
    """
    set -euo pipefail

    seqkit stats -j ${task.cpus} ${r1} ${r2} > ${sampleid}.${label}.stats.tsv
    """
}

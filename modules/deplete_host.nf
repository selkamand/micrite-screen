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

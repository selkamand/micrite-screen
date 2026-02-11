nextflow.enable.dsl = 2

// Required inputs
params.ref = null
// bowtie2 index prefix
params.kraken_db = null
// krakenuniq database directory
params.bam = null
// input bam
params.decoys = null
// decoy contigs list (required by your script)

// Minimal knobs (can be hardcoded if you want even fewer)
params.threads = 8
params.preload_size = '0'
params.outdir = 'results'

process FETCH_UNMAPPED_PRE {
    tag "${prefix}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(bam), path(decoys), val(prefix)

    output:
    tuple val(prefix), path("${prefix}_R1.fastq.gz"), path("${prefix}_R2.fastq.gz")

    script:
    """
  fetch_unmapped_reads.sh \
    --bam ${bam} \
    --decoys ${decoys} \
    --prefix ${prefix} \
    --threads ${params.threads} \
    --outdir .
  """
}

process FETCH_UNMAPPED_POST {
    tag "${prefix}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(bam), path(decoys), val(prefix)

    output:
    tuple val(prefix), path("${prefix}_R1.fastq.gz"), path("${prefix}_R2.fastq.gz")

    script:
    """
  fetch_unmapped_reads.sh \
    --bam ${bam} \
    --decoys ${decoys} \
    --prefix ${prefix} \
    --threads ${params.threads} \
    --outdir .
  """
}

process ALIGN_BOWTIE2 {
    tag "${prefix}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(prefix), path(r1), path(r2)

    output:
    tuple val(prefix), path("${prefix}.sorted.bam"), path("${prefix}.sorted.bam.bai")

    script:
    """
  align_bowtie2.sh \
    -x ${params.ref} \
    -1 ${r1} \
    -2 ${r2} \
    -o ${prefix} \
    -t ${params.threads} \
    -q 0
  """
}

process KRAKENUNIQ {
    tag "${prefix}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(prefix), path(r1), path(r2)

    output:
    path "${prefix}.krakenuniq.report.txt"

    script:
    """
  krakenuniq \
    --paired \
    --preload-size ${params.preload_size} \
    --threads ${params.threads} \
    --db ${params.kraken_db} \
    --report ${prefix}.krakenuniq.report.txt \
    --output off \
    ${r1} ${r2}
  """
}

workflow {

    // Fail fast: only required params
    if (!params.ref || !params.kraken_db || !params.bam || !params.decoys) {
        error(
            """
Missing required params:
  --ref <bowtie2 index prefix>
  --kraken_db <krakenuniq db dir>
  --bam <input.bam>
  --decoys <contigs.txt>
"""
        )
    }

    def bam = file(params.bam)
    def decoys = file(params.decoys)
    def sample = bam.baseName

    // 1) unmapped from original BAM
    unmapped1 = FETCH_UNMAPPED_PRE(channel.of(tuple(bam, decoys, "${sample}.pre")))

    // 2) align to reference
    aligned = ALIGN_BOWTIE2(unmapped1)

    // 3) unmapped from new BAM
    unmapped2 = aligned.map { prefix, bam2, bai2 -> tuple(bam2, decoys, "${sample}.post") }
        | FETCH_UNMAPPED_POST

    // 4) krakenuniq on final unmapped reads
    KRAKENUNIQ(unmapped2)
}

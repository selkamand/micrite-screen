nextflow.enable.dsl = 2

// [[ Required inputs ]]

// bowtie2 index prefix
params.ref = null

// krakenuniq database directory
params.kraken_db = null

// input bam
params.bam = null

// path to txt file with decoy contig names
params.decoys = null


// Optional
params.threads = 8
params.threads_kraken = 2
params.preload_size = '0'
params.outdir = 'results'

// Configure bowtie2 preset, default to sensitive 
// Allowed: very-fast, fast, sensitive, very-sensitive
params.bowtie2_preset = 'sensitive'

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
    --preset ${params.bowtie2_preset} \
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
    --threads ${params.threads_kraken} \
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

    // NEW: Validate preset early (so failures are obvious)
    def allowed_presets = ['very-fast', 'fast', 'sensitive', 'very-sensitive'] as Set
    if (!allowed_presets.contains(params.bowtie2_preset as String)) {
        error(
            """
Invalid --bowtie2_preset '${params.bowtie2_preset}'.
Valid values: very-fast, fast, sensitive, very-sensitive
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

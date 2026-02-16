nextflow.enable.dsl = 2

// Help flag (so `--help` works)
params.help = false

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
params.preload_size = '16G'
params.outdir = 'results'

// Configure bowtie2 preset, default to sensitive
// Allowed: very-fast, fast, sensitive, very-sensitive
params.bowtie2_preset = 'sensitive'


process FETCH_UNMAPPED_PRE {
    tag "${prefix}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(bam), path(bai), path(decoys), val(prefix)

    output:
    tuple val(prefix), path("${prefix}.R1.fq.gz"), path("${prefix}.R2.fq.gz")

    script:
    """
    ls;
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
    tuple val(prefix), path("${prefix}.R1.fq.gz"), path("${prefix}.R2.fq.gz")

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
    --preset ${params.bowtie2_preset}

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

    // Print help/usage and exit
    if (params.help) {
        log.info(
            """
USAGE:
  nextflow run main.nf \\
    --ref <bowtie2_index_prefix> \\
    --kraken_db <krakenuniq_db_dir> \\
    --bam <input.bam> \\
    --decoys <contigs.txt> \\
    [options]

REQUIRED:
  --ref            Bowtie2 index prefix for T2T reference genome to use in host depletion step (e.g. /path/to/index/prefix)
  --kraken_db      Krakenuniq database directory
  --bam            Input BAM file
  --decoys         Path to txt file with decoy contig names (one per line)

OPTIONAL (defaults shown):
  --outdir         Output directory (default: ${params.outdir})
  --threads        Threads for fetch_unmapped_reads + bowtie2 (default: ${params.threads})
  --threads_kraken  Threads for krakenuniq (default: ${params.threads_kraken})
  --preload_size   Krakenuniq preload size (default: ${params.preload_size})
  --bowtie2_preset Bowtie2 preset: very-fast|fast|sensitive|very-sensitive (default: ${params.bowtie2_preset})

HOW TO SET PARAMS:
  Use space or '=' forms:
    --threads 16
    --threads=16

EXAMPLE:
  nextflow run main.nf -resume \\
    --ref /refs/hg38/bowtie2/hg38 \\
    --kraken_db /db/krakenuniq \\
    --bam sample.bam \\
    --decoys decoys.txt \\
    --outdir results \\
    --threads 16 \\
    --threads_kraken 8 \\
    --preload_size 32G \\
    --bowtie2_preset very-sensitive
"""
        )
        System.exit(0)
    }

    // Fail fast: only required params
    if (!params.ref || !params.kraken_db || !params.bam || !params.decoys) {
        error(
            """
Missing required params:
  --ref <bowtie2 index prefix for T2T reference genome>
  --kraken_db <krakenuniq db dir>
  --bam <input.bam>
  --decoys <contigs.txt>

Tip: run with --help for full usage.
"""
        )
    }

    // Validate preset early (so failures are obvious)
    def allowed_presets = ['very-fast', 'fast', 'sensitive', 'very-sensitive'] as Set
    if (!allowed_presets.contains(params.bowtie2_preset as String)) {
        error(
            """
Invalid --bowtie2_preset '${params.bowtie2_preset}'.
Valid values: very-fast, fast, sensitive, very-sensitive

Tip: run with --help for full usage.
"""
        )
    }

    // Setup Paramaters
    def bam = file(params.bam)
    def decoys = file(params.decoys)
    def sample = bam.baseName

    // Check if required files exist
    if (!bam.exists()) {
        error("BAM file not found: ${bam}")
    }

    // define index path:
    def bai = file("${bam}.bai")

    // sample.bam -> sample.bam.bai
    if (!bai.exists()) {
        error("BAM index not found: ${bai}")
    }
    else {
        log.info("Bam index found: ${bai}")
    }

    if (!decoys.exists()) {
        error("Decoys file not found: ${decoys}")
    }


    ref = file(params.ref)
    bowtie_index = [
        file("${ref}.1.bt2"),
        file("${ref}.2.bt2"),
        file("${ref}.3.bt2"),
        file("${ref}.4.bt2"),
        file("${ref}.rev.1.bt2"),
        file("${ref}.rev.2.bt2"),
    ]


    if (!ref.exists()) {
        error("Reference genome not found ${ref}")
    }


    // Ensure krakenuniq database exists:
    if (!file(params.kraken_db).exists()) {
        error("Krakenuniq DB dir not found: ${params.kraken_db}")
    }

    // 1) unmapped from original BAM
    unmapped1 = FETCH_UNMAPPED_PRE(channel.of(tuple(bam, bai, decoys, "${sample}.pre")))

    // 2) align to reference
    aligned = ALIGN_BOWTIE2(unmapped1)

    // 3) unmapped from new BAM
    unmapped2 = aligned.map { prefix, bam2, bai2 -> tuple(bam2, decoys, "${sample}.post") }
        | FETCH_UNMAPPED_POST

    // 4) krakenuniq on final unmapped reads
    KRAKENUNIQ(unmapped2)
}

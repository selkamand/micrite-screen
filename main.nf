params {

    // Help flag (so `--help` works)
    help: Boolean = false

    // [[ Required inputs ]]
    // bowtie2 index prefix
    ref: Path

    // krakenuniq database directory
    kraken_db: Path

    // input bam
    bam: Path

    // path to txt file with decoy contig names
    decoys: Path

    // Sample Identified
    sampleid: String

    // Optional
    threads: Integer = 8
    threads_kraken: Integer = 2
    preload_size: String = "2G"

    outdir: Path = "micritescreen"

    // Configure bowtie2 preset, default to sensitive
    // Allowed: very-fast, fast, sensitive, very-sensitive
    bowtie2_preset: String = "sensitive"
}



process FETCH_UNMAPPED_PRE {
    input:
    tuple val(sampleid), path(bam), path(decoys), path(bai)

    output:
    tuple val(sampleid), path("${sampleid}.unmapped.R1.fq.gz"), path("${sampleid}.unmapped.R2.fq.gz")

    script:
    """
    ls;
     fetch_unmapped_reads.sh \
    --bam ${bam} \
    --decoys ${decoys} \
    --prefix "${sampleid}.unmapped" \
    --threads ${params.threads} \
    --outdir .
  """
}

process FETCH_UNMAPPED_POST {
    tag "Unmap from ${bam}"

    input:
    tuple val(sampleid), path(bam), path(decoys), path(bai)

    output:
    tuple val(sampleid), path("${sampleid}.hostdepleted.R1.fq.gz"), path("${sampleid}.hostdepleted.R2.fq.gz")

    script:
    """
    fetch_unmapped_reads.sh \
    --bam ${bam} \
    --decoys ${decoys} \
    --prefix ${sampleid}.hostdepleted \
    --threads ${params.threads} \
    --outdir .
  """
}

process ALIGN_BOWTIE2 {
    input:
    tuple val(sampleid), path(ref), path(r1), path(r2), path(bowtie_indexes)

    output:
    tuple val(sampleid), path("${sampleid}.realigned.sorted.bam"), path("${sampleid}.realigned.sorted.bam.bai")

    script:
    """
  align_bowtie2.sh \
    -x ${ref} \
    -1 ${r1} \
    -2 ${r2} \
    -o ${sampleid}.realigned \
    -t ${params.threads} \
    --preset ${params.bowtie2_preset}
  """
}

process KRAKENUNIQ {
    tag "${sampleid}"
    debug true

    input:
    tuple val(sampleid), path(krakendb), path(r1), path(r2)

    output:
    tuple val(sampleid), path("${sampleid}.krakenuniq.report.txt"), path("${sampleid}.kout.txt")

    script:
    """
    set -euo pipefail

    krakenuniq \
        --paired \
        --preload-size ${params.preload_size} \
        --threads ${params.threads_kraken} \
        --db ${krakendb} \
        --report ${sampleid}.krakenuniq.report.txt \
        --output ${sampleid}.kout.txt \
        ${r1} ${r2}
    """
}


process HOST_DEPLETION_STATS {
    input:
    tuple val(sampleid), path(bam_original), path(r1_unmapped), path(r2_unmapped), path(r1_depleted), path(r2_depleted), path(bai_original)

    output:
    tuple val(sampleid), path("${sampleid}.original.bam.stats.tsv"), path("${sampleid}.unmapped.stats.tsv"), path("${sampleid}.host_depleted.stats.tsv")

    script:
    """
    set -euo pipefail

    samtools flagstat ${bam_original} > ${sampleid}.original.bam.stats.tsv

    seqkit stats ${r1_unmapped} ${r2_unmapped} > ${sampleid}.unmapped.stats.tsv

    seqkit stats ${r1_depleted} ${r2_depleted} > ${sampleid}.host_depleted.stats.tsv
    """
}

workflow {

    main:
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
                --bowtie2_preset sensitive
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
    def kraken_db = file(params.kraken_db)
    def sampleid = params.sampleid

    // Check if required files exist
    if (!bam.exists()) {
        error("BAM file not found: ${bam}")
    }


    // define index paths:
    def bai = file("${bam}.bai")
    def ref = file(params.ref)

    def bowtie_index = [
        file("${ref}.1.bt2"),
        file("${ref}.2.bt2"),
        file("${ref}.3.bt2"),
        file("${ref}.4.bt2"),
        file("${ref}.rev.1.bt2"),
        file("${ref}.rev.2.bt2"),
    ]


    // File Checks

    // Ensure krakenuniq database exists:
    if (!file(params.kraken_db).exists()) {
        error("Krakenuniq DB dir not found: ${params.kraken_db}")
    }

    // Check Sample Identifier is Valid
    if (!(sampleid ==~ /^[A-Za-z0-9._-]+$/)) {
        error(
            """
            Invalid sampleid: '${params.sampleid}'

            Only letters, numbers, '.', '_', and '-' are allowed.
            No spaces or special characters.
            """
        )
    }

    // Check all required files exist
    bowtie_index.each { f ->
        if (!f.exists()) {
            error("Missing bowtie index file: ${f}")
        }
    }

    if (!bai.exists()) {
        error("Missing required bam index file: ${bai}")
    }

    if (!ref.exists()) {
        error("Failed to find reference genome for host depletion at: ${ref}")
    }


    // 1) unmapped from original BAM
    unmapped_reads1_ch = FETCH_UNMAPPED_PRE(channel.of(tuple(sampleid, bam, decoys, bai)))


    // 2) align to reference
    realigned_ch = unmapped_reads1_ch.map { s, r1, r2 -> tuple(s, ref, r1, r2, bowtie_index) }
        | ALIGN_BOWTIE2


    // 3) unmapped from new BAM
    host_depleted_ch = realigned_ch.map { s, bam_realigned, bai_realigned -> tuple(s, bam_realigned, decoys, bai_realigned) }
        | FETCH_UNMAPPED_POST

    // 4) krakenuniq on final unmapped reads
    host_depleted_ch.map { s, r1, r2 -> tuple(s, kraken_db, r1, r2) }
        | KRAKENUNIQ

    // 5) stats
    // Channel carrying the "original" context we need for stats
    orig_ch = channel.of(tuple(sampleid, bam, bai))

    // Join everything by sampleid (tuple index 0), then reorder into HOST_DEPLETION_STATS input format
    stats_in_ch = orig_ch
        .join(unmapped_reads1_ch, by: 0)
        .join(host_depleted_ch, by: 0)
        .map { sid, bam0, bai0, r1u, r2u, r1d, r2d ->
            tuple(sid, bam0, r1u, r2u, r1d, r2d, bai0)
        }

    // Now the process will actually run
    host_depletion_stats_ch = HOST_DEPLETION_STATS(stats_in_ch)

    publish:
    host_depleted_reads = FETCH_UNMAPPED_POST.out
    krakenuniq = KRAKENUNIQ.out
    stats = HOST_DEPLETION_STATS.out
}

// Outputs to save in final directory
output {
    host_depleted_reads {
        path "${params.outdir}/${params.sampleid}"
        mode 'copy'
    }
    krakenuniq {
        path "${params.outdir}/${params.sampleid}"
        mode 'copy'
    }
    stats {
        path "${params.outdir}/${params.sampleid}stats"
        mode 'copy'
    }
}

#!/usr/bin/env nextflow

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
    outdir: Path = "micritescreen"

    // Configure bowtie2 preset, default to sensitive
    // Allowed: very-fast, fast, sensitive, very-sensitive
    bowtie2_preset: String = "sensitive"
}

include { DEPLETE_HOST } from "./subworkflows/host_depletion.nf"
include { KRAKENUNIQ } from './modules/classify.nf'

def usage() {
    log.info(
        """
    USAGE:
      nextflow run main.nf \\
        --ref <bowtie2_index_prefix> \\
        --bam <input.bam> \\
        --decoys <contigs.txt> \\
        --sampleid <sample_id> \\
        [--kraken_db <krakenuniq_db_dir>] \\
        [options]

    REQUIRED:
      --ref              Bowtie2 index prefix
      --bam              Input BAM
      --decoys           Decoy contig list
      --sampleid         Sample identifier

    OPTIONAL:
      --kraken_db            Krakenuniq DB directory; required only if --run_classification true
      --outdir               Output directory (default: ${params.outdir})
      --bowtie2_preset       very-fast|fast|sensitive|very-sensitive (default: ${params.bowtie2_preset})
    """
    )
    System.exit(0)
}

def validateParams() {
    def missing = []

    if (!params.ref) {
        missing << '--ref'
    }
    if (!params.bam) {
        missing << '--bam'
    }
    if (!params.decoys) {
        missing << '--decoys'
    }
    if (!params.sampleid) {
        missing << '--sampleid'
    }

    if (missing) {
        error("Missing required params: ${missing.join(', ')}\nTip: run with --help")
    }

    def allowed_presets = ['very-fast', 'fast', 'sensitive', 'very-sensitive'] as Set
    if (!allowed_presets.contains(params.bowtie2_preset as String)) {
        error("Invalid --bowtie2_preset '${params.bowtie2_preset}'. Valid values: ${allowed_presets.join(', ')}")
    }

    if (!(params.sampleid ==~ /^[A-Za-z0-9._-]+$/)) {
        error("Invalid sampleid '${params.sampleid}'. Only letters, numbers, '.', '_', and '-' are allowed.")
    }
}

def validateFiles() {
    def bam = file(params.bam)
    def bai = file("${bam}.bai")
    def decoys = file(params.decoys)
    def ref = file(params.ref)

    if (!bam.exists()) {
        error("BAM file not found: ${bam}")
    }
    if (!bai.exists()) {
        error("Missing BAM index: ${bai}")
    }
    if (!decoys.exists()) {
        error("Decoy list not found: ${decoys}")
    }
    if (!ref.exists()) {
        error("Reference prefix not found: ${ref}")
    }

    def bowtie_index = [
        file("${ref}.1.bt2"),
        file("${ref}.2.bt2"),
        file("${ref}.3.bt2"),
        file("${ref}.4.bt2"),
        file("${ref}.rev.1.bt2"),
        file("${ref}.rev.2.bt2"),
    ]

    bowtie_index.each { f ->
        if (!f.exists()) {
            error("Missing bowtie2 index file: ${f}")
        }
    }

    if (params.run_classification) {
        def kraken_db = file(params.kraken_db)
        if (!kraken_db.exists()) {
            error("Krakenuniq DB dir not found: ${kraken_db}")
        }
    }
}

def buildBowtieIndex(ref) {
    [
        file("${ref}.1.bt2"),
        file("${ref}.2.bt2"),
        file("${ref}.3.bt2"),
        file("${ref}.4.bt2"),
        file("${ref}.rev.1.bt2"),
        file("${ref}.rev.2.bt2"),
    ]
}


workflow {

    main:
    // Print help/usage and exit
    if (params.help) {
        usage()
    }

    // Fail fast if params are weird or files need validating
    validateParams()
    validateFiles()

    // Setup Paramaters
    def bam = file(params.bam)
    def bai = file("${bam}.bai")
    def decoys = file(params.decoys)
    def kraken_db = file(params.kraken_db)
    def ref = file(params.ref)
    def bowtie_index = [
        file("${ref}.1.bt2"),
        file("${ref}.2.bt2"),
        file("${ref}.3.bt2"),
        file("${ref}.4.bt2"),
        file("${ref}.rev.1.bt2"),
        file("${ref}.rev.2.bt2"),
    ]

    host_depletion_input_ch = channel.of(
        tuple(
            params.sampleid,
            bam,
            bai,
            decoys,
            ref,
            bowtie_index,
            params.bowtie2_preset,
        )
    )


    // Deplete Host Reads
    host_depleted_ch = DEPLETE_HOST(host_depletion_input_ch)

    // Run KrakenUniq to classify reads 
    krakenuniq_ch = host_depleted_ch.reads.map { sid, r1, r2 ->
        tuple(sid, kraken_db, r1, r2)
    }
        | KRAKENUNIQ

    publish:
    host_depleted_reads = host_depleted_ch.reads
    stats = host_depleted_ch.stats
    krakenuniq = krakenuniq_ch
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
        path "${params.outdir}/${params.sampleid}/stats"
        mode 'copy'
    }
}

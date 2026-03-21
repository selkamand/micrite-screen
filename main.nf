#!/usr/bin/env nextflow

params {

    // Help flag (so `--help` works)
    help: Boolean = false

    // Run host depletion and microbial screen on 'bam' or 'fastq'
    mode = "bam"

    // Input bam (reads aligned to host genome)
    // Either a bam or a pair of fastq files must be supplied
    bam: Path?

    // Input FASTQs (Must be supplied if mode = "fastq") 
    r1: Path?
    r2: Path?

    // [[ Required inputs ]]
    // bowtie2 index prefix
    ref: Path

    // krakenuniq database directory
    kraken_db: Path

    // Path to txt file with decoy contig names
    decoys: Path

    // Sample Identified
    sampleid: String

    // Optional
    outdir: Path = "micritescreen"

    // Configure bowtie2 preset, default to sensitive
    // Allowed: very-fast, fast, sensitive, very-sensitive
    bowtie2_preset: String = "sensitive"

    // Run pipeline up to host-depleted reads but do not proceed to classification
    skip_classification: Boolean = false
}

include { DEPLETE_HOST_FASTQ } from "./subworkflows/host_depletion.nf"
include { FETCH_UNMAPPED } from "./modules/deplete_host.nf"
include { KRAKENUNIQ } from './modules/classify.nf'
include { KREPORT_TO_KRONA } from "./modules/kreport_to_krona.nf"
include { BAMSTATS } from "./modules/stats.nf"

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
      --kraken_db            Krakenuniq DB directory; required only if --skip_classification is false
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
    if (!params.bam & params.mode == "bam") {
        missing << '--bam'
    }
    if (!params.decoys) {
        missing << '--decoys'
    }
    if (!params.sampleid) {
        missing << '--sampleid'
    }
    if (!params.r1 & params.mode == "fastq") {
        missing << '--r1'
    }
    if (!params.r2 & params.mode == "fastq") {
        missing << '--r2'
    }

    if (missing) {
        error("Missing required params: ${missing.join(', ')}\nTip: run with --help")
    }

    def allowed_modes = ['bam', 'fastq'] as Set

    if (!allowed_modes.contains(params.mode as String)) {
        error("Invalid --mode '${params.mode}'. Valid values: ${allowed_modes.join(', ')}")
    }

    def allowed_presets = ['very-fast', 'fast', 'sensitive', 'very-sensitive'] as Set

    if (!allowed_presets.contains(params.bowtie2_preset as String)) {
        error("Invalid --bowtie2_preset '${params.bowtie2_preset}'. Valid values: ${allowed_presets.join(', ')}")
    }

    if (!(params.sampleid ==~ /^[A-Za-z0-9._-]+$/)) {
        error("Invalid sampleid '${params.sampleid}'. Only letters, numbers, '.', '_', and '-' are allowed.")
    }

    if (params.mode == "bam") {
        if (params.r1) {
            error("There is no reason to supply fastqs to --r1 when --mode = 'bam'. Drop --r1 parameter or change --mode to 'fastq'")
        }
        if (params.r2) {
            error("There is no reason to supply fastqs to --r2 when --mode = 'bam'. Drop --r2 parameter or change --mode to 'fastq'")
        }
    }
}

def validateFiles() {

    def mode = params.mode

    if (mode == "bam") {
        def bam = file(params.bam)
        def bai = file("${bam}.bai")
        if (!bam.exists()) {
            error("BAM file not found: ${bam}")
        }
        if (!bai.exists()) {
            error("Missing BAM index: ${bai}")
        }
    }
    if (mode == "fastq") {
        def r1 = file(params.r1)
        def r2 = file(params.r2)


        if (!r1.exists()) {
            error("FASTQ file not found: ${r1}")
        }
        if (!r2.exists()) {
            error("FASTQ file not found: ${r2}")
        }
    }

    def decoys = file(params.decoys)
    def ref = file(params.ref)

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

    if (!params.skip_classification) {
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

    // Setup Reference Genome Paramaters
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

    // Setup refindex channel
    ch_host_refgenome = channel.of(tuple(ref, bowtie_index, params.bowtie2_preset))
    ch_decoys = channel.of(decoys)


    ch_input_fastqs = channel.empty()
    ch_unmapped_fastqs = channel.empty()
    ch_original_bamstats = channel.empty()

    // To run host depletion from BAM we first need to extract unmapped reads into ch_input_fastqs
    if (params.mode == "bam") {
        def bam = file(params.bam)
        def bai = file("${bam}.bai")

        // Add stats of bam. Rename st 
        ch_input_bam = channel.of(tuple(params.sampleid, bam, bai))
        ch_original_bamstats = BAMSTATS("original", ch_input_bam)
        ch_input_fastqs = FETCH_UNMAPPED("unmapped", ch_input_bam, ch_decoys)
        ch_unmapped_fastqs = ch_input_fastqs
    }
    else if (params.mode == "fastq") {
        def r1 = file(params.r1)
        def r2 = file(params.r2)

        ch_input_fastqs = channel.of(
            tuple(
                params.sampleid,
                r1,
                r2,
            )
        )
    }
    else {
        error("micrite-screen does not currently support param.mode = [${params.mode}]")
    }

    // Deplete host from reads
    ch_host_depleted = DEPLETE_HOST_FASTQ(ch_input_fastqs, ch_host_refgenome, ch_decoys)

    // Always define kraken channels so publish: can see them even if classification is skipped
    krakenuniq_ch = channel.empty()
    krona_ch = channel.empty()

    if (!params.skip_classification) {
        // Run KrakenUniq to classify reads 
        krakenuniq_ch = ch_host_depleted.reads.map { sid, read1, read2 ->
            tuple(sid, kraken_db, read1, read2)
        }
            | KRAKENUNIQ

        // Generate Krona compatible file from report
        krona_ch = KREPORT_TO_KRONA(krakenuniq_ch)
    }

    publish:
    original_bamstats = ch_original_bamstats
    unmapped_reads = ch_unmapped_fastqs
    host_depleted_reads = ch_host_depleted.reads
    depletion_stats = ch_host_depleted.stats
    krakenuniq = krakenuniq_ch
    krona = krona_ch
}

// Outputs to save in final directory
output {
    original_bamstats {
        path "${params.outdir}/${params.sampleid}/stats"
        mode 'copy'
    }
    unmapped_reads {
        path "${params.outdir}/${params.sampleid}/reads"
        mode 'copy'
    }
    host_depleted_reads {
        path "${params.outdir}/${params.sampleid}/reads"
        mode 'copy'
    }
    depletion_stats {
        path "${params.outdir}/${params.sampleid}/stats"
        mode 'copy'
    }
    krakenuniq {
        path "${params.outdir}/${params.sampleid}"
        mode 'copy'
    }
    krona {
        path "${params.outdir}/${params.sampleid}/krona"
        mode 'copy'
    }
}

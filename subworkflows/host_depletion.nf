include { FETCH_UNMAPPED ; ALIGN_BOWTIE2 ; HOST_DEPLETION_STATS_FASTQ } from '../modules/deplete_host.nf'

workflow DEPLETE_HOST_FASTQ {
    take:
    ch_reads // channel: [ val(sampleid), path(r1), path(r2) ]
    ch_reference_genome // channel: [ path(reference_genome), [path(bowtie_index1),path(bowtie_index2), ...]],
    ch_decoys // channel: path(decoy_names)

    main:

    // Realign unmapped reads to T2T genome
    ch_realigned_bam = ALIGN_BOWTIE2("realigned", ch_reads, ch_reference_genome)

    // Fetch Unmapped reads from T2T realignment
    ch_host_depleted_reads = FETCH_UNMAPPED("host_depleted", ch_realigned_bam, ch_decoys)

    // Compute Stats
    ch_stats_in = ch_reads.join(ch_host_depleted_reads, by: 0)
    stats_ch = HOST_DEPLETION_STATS_FASTQ(ch_stats_in)

    emit:
    reads = ch_host_depleted_reads //tuple val(sampleid), path("${sampleid}.hostdepleted.R1.fq.gz"), path("${sampleid}.hostdepleted.R2.fq.gz")
    stats = stats_ch // tuple val(sampleid), path("${sampleid}.original.stats.tsv"), path("${sampleid}.host_depleted.stats.tsv")
}

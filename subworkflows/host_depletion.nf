include { FETCH_UNMAPPED ; ALIGN_BOWTIE2 } from '../modules/deplete_host.nf'
include { BAMSTATS ; SEQKIT_STATS as SEQKIT_STATS_PREALIGNMENT ; SEQKIT_STATS as SEQKIT_STATS_POST_DEPLETION } from '../modules/stats.nf'


workflow DEPLETE_HOST_FASTQ {
    take:
    ch_reads // channel: [ val(sampleid), path(r1), path(r2) ]
    ch_reference_genome // channel: [ path(reference_genome), [path(bowtie_index1),path(bowtie_index2), ...]],
    ch_decoys // channel: path(decoy_names)

    main:

    // Realign unmapped reads to T2T genome
    ch_realigned_bam = ALIGN_BOWTIE2("realigned", ch_reads, ch_reference_genome)
    ch_realigned_bamstats = BAMSTATS("realigned", ch_realigned_bam)

    // Fetch Unmapped reads from T2T realignment
    ch_host_depleted_reads = FETCH_UNMAPPED("host_depleted", ch_realigned_bam, ch_decoys)

    // Compute Read Stats
    ch_stats_before_realignment = SEQKIT_STATS_PREALIGNMENT("pre_realignment", ch_reads)
    ch_stats_after_host_depleted = SEQKIT_STATS_POST_DEPLETION("host_depleted", ch_host_depleted_reads)
    ch_stats = ch_stats_before_realignment.join(ch_stats_after_host_depleted).join(ch_realigned_bamstats)

    emit:
    reads = ch_host_depleted_reads //tuple val(sampleid), path("${sampleid}.hostdepleted.R1.fq.gz"), path("${sampleid}.hostdepleted.R2.fq.gz")
    stats = ch_stats // tuple val(sampleid), path("${sampleid}.pre_realignment.stats.tsv"), path("${sampleid}.host_depleted.stats.tsv")
}

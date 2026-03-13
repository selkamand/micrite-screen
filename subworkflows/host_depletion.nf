include { FETCH_UNMAPPED_PRE ; FETCH_UNMAPPED_POST ; FETCH_UNMAPPED ; ALIGN_BOWTIE2 ; ALIGN_BOWTIE2_NEW ; HOST_DEPLETION_STATS_BAM ; HOST_DEPLETION_STATS_FASTQ } from '../modules/deplete_host.nf'

workflow DEPLETE_HOST_BAM {
    take:
    host_depletion_input_ch

    main:

    // Original Bam Information
    original_bam_ch = host_depletion_input_ch.map { sid, bam, bai, decoys, _ref, _bowtie_index, _preset ->
        tuple(sid, bam, bai, decoys)
    }

    // Unmap Reads
    unmapped_reads_ch = FETCH_UNMAPPED_PRE(original_bam_ch)

    // Prep Bowtie Input
    bowtie_input_ch = unmapped_reads_ch
        .join(host_depletion_input_ch, by: 0)
        .map { sid, r1, r2, _bam, _bai, _decoys, ref, bowtie_index, preset ->
            tuple(sid, ref, r1, r2, bowtie_index, preset)
        }

    // Realign unmapped reads to T2T genome
    realigned_bam_ch = ALIGN_BOWTIE2(bowtie_input_ch)

    // Prep input for final host depletion (unmapping from T2T) 
    second_unmap_input_ch = realigned_bam_ch
        .join(host_depletion_input_ch, by: 0)
        .map { sid, bam_realigned, bai_realigned, _bam_original, _bai_original, decoys, _ref, _bowtie_index, _preset ->
            tuple(
                sid,
                bam_realigned,
                bai_realigned,
                decoys,
            )
        }

    // Fetch Unmapped reads from T2T realignment
    host_depleted_reads_ch = FETCH_UNMAPPED_POST(second_unmap_input_ch)

    // Compute Stats
    stats_in_ch = original_bam_ch
        .join(unmapped_reads_ch, by: 0)
        .join(host_depleted_reads_ch, by: 0)

    stats_ch = HOST_DEPLETION_STATS_BAM(stats_in_ch)

    emit:
    reads = host_depleted_reads_ch
    stats = stats_ch
}

workflow DEPLETE_HOST_FASTQ {
    take:
    ch_reads // channel: [ val(sampleid), path(r1), path(r2) ]
    ch_reference_genome // channel: [ path(reference_genome), [path(bowtie_index1),path(bowtie_index2), ...]],
    ch_decoys // channel: path(decoy_names)

    main:

    // Prep Bowtie Input
    //bowtie_input_ch = fastq_to_deplete_ch

    // Realign unmapped reads to T2T genome
    realigned_bam_ch = ALIGN_BOWTIE2_NEW(ch_reads, ch_reference_genome)

    // Prep input for final host depletion (unmapping from T2T) 
    second_unmap_input_ch = realigned_bam_ch.map { sid, bam_realigned, bai_realigned, decoys, _ref, _bowtie_index, _preset ->
        tuple(
            sid,
            bam_realigned,
            bai_realigned,
            decoys,
        )
    }

    // Fetch Unmapped reads from T2T realignment
    ch_host_depleted_reads = FETCH_UNMAPPED(second_unmap_input_ch, ch_decoys)

    // Compute Stats
    ch_stats_in = ch_reads.join(ch_host_depleted_reads, by: 0)

    stats_ch = HOST_DEPLETION_STATS_FASTQ(ch_stats_in)

    emit:
    reads = ch_host_depleted_reads
    stats = stats_ch
}

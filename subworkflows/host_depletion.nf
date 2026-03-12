include { FETCH_UNMAPPED_PRE ; FETCH_UNMAPPED_POST ; ALIGN_BOWTIE2 ; HOST_DEPLETION_STATS } from '../modules/deplete_host.nf'

workflow DEPLETE_HOST {
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
        .view()
        .map { sid, r1, r2, _bam, _bai, _decoys, ref, bowtie_index, preset ->
            tuple(sid, ref, r1, r2, bowtie_index, preset)
        }

    // Realign unmapped reads to T2T genome
    realigned_bam_ch = ALIGN_BOWTIE2(bowtie_input_ch)

    // Prep input for final host depletion (unmapping from T2T) 
    second_unmap_input_ch = realigned_bam_ch
        .join(host_depletion_input_ch, by: 0)
        .map { sid, bam_realigned, bai_realigned, decoys, _ref, _bowtie_index, _preset ->
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
        .join(original_bam_ch, by: 0)
        .join(host_depleted_reads_ch, by: 0)

    stats_ch = HOST_DEPLETION_STATS(stats_in_ch)

    emit:
    reads = host_depleted_reads_ch
    stats = stats_ch
}

process KRAKENUNIQ {
    tag "${sampleid}"
    debug true
    cpus params.threads_kraken
    ext preload_size: params.preload_size

    input:
    tuple val(sampleid), path(krakendb), path(r1), path(r2)

    output:
    tuple val(sampleid), path("${sampleid}.krakenuniq.report.txt"), path("${sampleid}.kout.txt")

    script:
    """
    set -euo pipefail

    krakenuniq \
        --paired \
        --preload-size ${task.ext.preload_size} \
        --threads ${task.cpus} \
        --db ${krakendb} \
        --report ${sampleid}.krakenuniq.report.txt \
        --output ${sampleid}.kout.txt \
        ${r1} ${r2}
    """
}

nextflow.enable.dsl = 2

process KRAKENUNIQ {
    tag "${sampleid}"
    debug true

    input:
    tuple val(sampleid), path(krakendb), path(r1), path(r2)

    output:
    tuple val(sampleid), path("${sampleid}.krakenuniq.report.txt"), path("${sampleid}.kout.txt")

    script:
    def preload_size = task.ext.args ?: '2G'
    """
    set -euo pipefail

    krakenuniq \
        --paired \
        --preload-size ${preload_size} \
        --threads ${task.cpus} \
        --db ${krakendb} \
        --report ${sampleid}.krakenuniq.report.txt \
        --output ${sampleid}.kout.txt \
        ${r1} ${r2}
    """
}

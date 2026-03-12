nextflow.enable.dsl = 2

process KREPORT_TO_KRONA {

    container "selkamandcci/micrite-sleuth:0.0.1"
    tag "${sampleid}"
    debug true

    input:
    tuple val(sampleid), path(kreport), path(kout)

    output:
    tuple val(sampleid), path("report.nointermediate.krona"), path("report.krona")

    script:
    """
    set -euo pipefail
    kreport2krona.py --no-intermediate-ranks -r ${kreport} -o report.nointermediate.krona
    kreport2krona.py -r ${kreport} -o report.krona
    # ktUpdateTaxonomy.sh
    # ktImportTaxonomy -i -o report.nointermediate.krona.html report.nointermediate.krona
    # ktImportTaxonomy -i -o report.krona.html report.krona
    """
}

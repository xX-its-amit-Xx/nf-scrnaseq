process REPORT {

    tag "${meta.id}"

    input:
    tuple val(meta), path(clustered_h5ad), path(qc_json)

    output:
    tuple val(meta), path("${meta.id}.report.html"), emit: html

    script:
    """
    make_report.py \\
        --h5ad ${clustered_h5ad} \\
        --qc-json ${qc_json} \\
        --sample-id ${meta.id} \\
        --tissue ${meta.tissue} \\
        --template ${projectDir}/assets/report_template.html \\
        --out ${meta.id}.report.html
    """

    stub:
    """
    touch ${meta.id}.report.html
    """
}

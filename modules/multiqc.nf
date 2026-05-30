process MULTIQC {

    tag "aggregate"

    publishDir "${params.outdir}/multiqc", mode: params.publish_mode

    input:
    path  '*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data

    script:
    def config_arg = file(params.multiqc_config).exists() ? "--config ${params.multiqc_config}" : ''
    """
    multiqc ${config_arg} --force .
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_data
    """
}

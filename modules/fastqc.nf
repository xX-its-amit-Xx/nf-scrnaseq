process FASTQC {

    tag "${meta.id}_${meta.read}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${reads}
    """

    stub:
    """
    touch ${meta.id}_${meta.read}_fastqc.html
    touch ${meta.id}_${meta.read}_fastqc.zip
    """
}

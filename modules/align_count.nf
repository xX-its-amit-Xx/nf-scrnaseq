process ALIGN_COUNT {

    tag "${meta.id}"

    input:
    tuple val(meta), path(r1), path(r2)
    path  kb_index
    path  t2g

    output:
    tuple val(meta), path("${meta.id}/counts_unfiltered"), emit: counts
    tuple val(meta), path("${meta.id}/run_info.json"),     emit: logs
    tuple val(meta), path("${meta.id}/inspect.json"),      emit: inspect, optional: true

    script:
    """
    mkdir -p ${meta.id}

    kb count \\
        -i ${kb_index} \\
        -g ${t2g} \\
        -x ${params.chemistry} \\
        -o ${meta.id} \\
        --h5ad \\
        -t ${task.cpus} \\
        --overwrite \\
        ${r1} ${r2}
    """

    stub:
    """
    mkdir -p ${meta.id}/counts_unfiltered
    touch   ${meta.id}/counts_unfiltered/adata.h5ad
    echo '{"n_reads": 0}' > ${meta.id}/run_info.json
    """
}

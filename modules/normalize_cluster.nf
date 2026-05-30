process NORMALIZE_CLUSTER {

    tag "${meta.id}"

    input:
    tuple val(meta), path(qc_h5ad)

    output:
    tuple val(meta), path("${meta.id}.clustered.h5ad"),   emit: h5ad
    tuple val(meta), path("${meta.id}.cluster_plots/*"),  emit: plots, optional: true

    script:
    """
    normalize_cluster.py \\
        --in-h5ad ${qc_h5ad} \\
        --sample-id ${meta.id} \\
        --n-top-genes ${params.n_top_genes} \\
        --n-pcs ${params.n_pcs} \\
        --n-neighbors ${params.n_neighbors} \\
        --resolution ${params.resolution} \\
        --out-h5ad ${meta.id}.clustered.h5ad \\
        --out-plotdir ${meta.id}.cluster_plots
    """

    stub:
    """
    touch ${meta.id}.clustered.h5ad
    mkdir -p ${meta.id}.cluster_plots
    """
}

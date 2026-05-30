process SCANPY_QC {

    tag "${meta.id}"

    input:
    tuple val(meta), path(counts_dir)

    output:
    tuple val(meta), path("${meta.id}.qc.h5ad"),     emit: h5ad
    tuple val(meta), path("${meta.id}.qc.json"),     emit: qc_metrics
    tuple val(meta), path("${meta.id}.qc_plots/*"),  emit: plots, optional: true

    script:
    def mt_arg = params.mt_gene_ids ? "--mt-gene-ids ${params.mt_gene_ids}" : ''
    """
    scanpy_qc.py \\
        --counts ${counts_dir} \\
        --sample-id ${meta.id} \\
        --tissue ${meta.tissue} \\
        --min-genes ${params.min_genes} \\
        --min-cells ${params.min_cells} \\
        --max-mito-pct ${params.max_mito_pct} \\
        --doublet-thresh ${params.doublet_thresh} \\
        ${mt_arg} \\
        --out-h5ad ${meta.id}.qc.h5ad \\
        --out-json ${meta.id}.qc.json \\
        --out-plotdir ${meta.id}.qc_plots
    """

    stub:
    """
    touch ${meta.id}.qc.h5ad
    echo '{"n_cells_pre": 0, "n_cells_post": 0}' > ${meta.id}.qc.json
    mkdir -p ${meta.id}.qc_plots
    """
}

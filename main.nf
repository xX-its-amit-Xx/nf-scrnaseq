#!/usr/bin/env nextflow
/*
 * nf-scrnaseq: a Nextflow DSL2 pipeline for single-cell RNA-seq
 * From raw FASTQs to a QC'd, clustered AnnData object + HTML report.
 */

nextflow.enable.dsl = 2

include { FASTQC }             from './modules/fastqc.nf'
include { ALIGN_COUNT }        from './modules/align_count.nf'
include { SCANPY_QC }          from './modules/scanpy_qc.nf'
include { NORMALIZE_CLUSTER }  from './modules/normalize_cluster.nf'
include { REPORT }             from './modules/report.nf'
include { MULTIQC }            from './modules/multiqc.nf'

def helpMessage() {
    log.info """
    ===========================================================
     nf-scrnaseq  v${workflow.manifest.version ?: '0.1.0'}
    ===========================================================
    Usage:
      nextflow run . --samplesheet samplesheet.csv -profile docker

    Required:
      --samplesheet      CSV with columns: sample_id,fastq_1,fastq_2,tissue
      --kb_index         Pre-built kb-python index (.idx)
      --t2g              transcript-to-gene mapping (t2g.txt)

    Optional QC params:
      --min_genes        Minimum genes per cell (default: ${params.min_genes})
      --max_mito_pct     Maximum mitochondrial percent (default: ${params.max_mito_pct})
      --n_top_genes      Highly variable genes count (default: ${params.n_top_genes})
      --resolution       Leiden clustering resolution (default: ${params.resolution})

    Output:
      --outdir           Results directory (default: ${params.outdir})

    Profiles:
      -profile standard  Local execution, conda/native binaries
      -profile docker    Local execution, containerized
      -profile slurm     SLURM cluster execution
      -profile test      Tiny bundled dataset for CI / smoke test
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

workflow {

    if (!params.samplesheet) {
        error "Missing --samplesheet. See `nextflow run . --help`."
    }

    // ----------------------------------------------------------
    //  Parse samplesheet -> per-sample channel of [meta, R1, R2]
    // ----------------------------------------------------------
    ch_samples = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id:     row.sample_id,
                tissue: row.tissue ?: 'unspecified'
            ]
            def r1 = file(row.fastq_1, checkIfExists: true)
            def r2 = file(row.fastq_2, checkIfExists: true)
            tuple(meta, r1, r2)
        }

    // ----------------------------------------------------------
    //  Stage 1: FastQC on raw reads (per FASTQ)
    // ----------------------------------------------------------
    ch_fastqc_in = ch_samples
        .flatMap { meta, r1, r2 -> [
            tuple(meta + [read: 'R1'], r1),
            tuple(meta + [read: 'R2'], r2)
        ]}

    FASTQC( ch_fastqc_in )

    // ----------------------------------------------------------
    //  Stage 2: Alignment + count matrix (kb-python / kallisto-bustools)
    // ----------------------------------------------------------
    ch_index = file(params.kb_index, checkIfExists: true)
    ch_t2g   = file(params.t2g,      checkIfExists: true)

    ALIGN_COUNT(
        ch_samples,
        ch_index,
        ch_t2g
    )

    // ----------------------------------------------------------
    //  Stage 3: Scanpy QC + doublet detection
    // ----------------------------------------------------------
    SCANPY_QC(
        ALIGN_COUNT.out.counts
    )

    // ----------------------------------------------------------
    //  Stage 4: Normalize, HVG, PCA, neighbors, Leiden, UMAP
    // ----------------------------------------------------------
    NORMALIZE_CLUSTER(
        SCANPY_QC.out.h5ad
    )

    // ----------------------------------------------------------
    //  Stage 5: HTML report per sample
    // ----------------------------------------------------------
    REPORT(
        NORMALIZE_CLUSTER.out.h5ad
            .join(SCANPY_QC.out.qc_metrics)
    )

    // ----------------------------------------------------------
    //  Stage 6: MultiQC aggregate report
    // ----------------------------------------------------------
    ch_multiqc = FASTQC.out.zip.map { meta, z -> z }
        .mix( ALIGN_COUNT.out.logs.map { meta, l -> l } )
        .collect()

    MULTIQC( ch_multiqc )
}

workflow.onComplete {
    log.info """
    -----------------------------------------------------------
     Pipeline completed: ${workflow.success ? 'SUCCESS' : 'FAILED'}
     Duration  : ${workflow.duration}
     Results   : ${params.outdir}
    -----------------------------------------------------------
    """.stripIndent()
}

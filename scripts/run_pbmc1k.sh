#!/usr/bin/env bash
# Real-world recipe: 10x Genomics pbmc_1k_v3 -> nf-scrnaseq end-to-end.
#
# Uses kb-python's pre-built human reference to skip the 25-min index build.
# Total wall time on a 4-core / 16 GB Codespace: ~25-35 min.
#   - download FASTQs                  ~3 min  (~1.5 GB)
#   - download pre-built kb-python ref ~5 min  (~3 GB)
#   - alignment + count                ~15 min
#   - QC + cluster + report            ~5 min
#
# Idempotent: re-running this script skips already-downloaded artifacts and
# leverages Nextflow's -resume to skip already-computed pipeline stages.

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

DATA_DIR="${DATA_DIR:-data/pbmc1k}"
REF_DIR="${REF_DIR:-refs/human}"
SAMPLES_CSV="${SAMPLES_CSV:-pbmc1k.csv}"
IMAGE="${IMAGE:-nf-scrnaseq:0.1.0}"

mkdir -p "$DATA_DIR" "$REF_DIR"

# ----------------------------------------------------------
#  1. FASTQs — 10x Genomics public PBMC 1k (v3 chemistry)
# ----------------------------------------------------------
FASTQ_TAR="$DATA_DIR/pbmc_1k_v3_fastqs.tar"
FASTQ_R1="$DATA_DIR/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz"
FASTQ_R2="$DATA_DIR/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz"

if [[ ! -f "$FASTQ_R1" || ! -f "$FASTQ_R2" ]]; then
    echo ">>> Downloading 10x pbmc_1k_v3 FASTQs (~1.5 GB)"
    curl -L --fail --retry 3 -o "$FASTQ_TAR" \
        "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar"
    tar xf "$FASTQ_TAR" -C "$DATA_DIR"
    rm -f "$FASTQ_TAR"
fi

# ----------------------------------------------------------
#  2. Reference — kb-python's pre-built human transcriptome index
#     `kb ref -d human` fetches a ready-made kallisto index + t2g.
# ----------------------------------------------------------
if [[ ! -f "$REF_DIR/index.idx" || ! -f "$REF_DIR/t2g.txt" ]]; then
    echo ">>> Downloading pre-built kb-python human reference (~3 GB)"
    docker run --rm -v "$PWD/$REF_DIR:/ref" -w /ref "$IMAGE" \
        kb ref -d human -i index.idx -g t2g.txt -f1 cdna.fa
fi

# ----------------------------------------------------------
#  3. Samplesheet
# ----------------------------------------------------------
cat > "$SAMPLES_CSV" <<EOF
sample_id,fastq_1,fastq_2,tissue
pbmc1k,$PWD/$FASTQ_R1,$PWD/$FASTQ_R2,pbmc
EOF
echo ">>> Samplesheet written to $SAMPLES_CSV"

# ----------------------------------------------------------
#  4. Run the pipeline
# ----------------------------------------------------------
echo ">>> Launching nf-scrnaseq"
nextflow run . \
    --samplesheet  "$SAMPLES_CSV" \
    --kb_index     "$PWD/$REF_DIR/index.idx" \
    --t2g          "$PWD/$REF_DIR/t2g.txt" \
    --chemistry    10xv3 \
    --outdir       results_pbmc1k \
    -profile       docker \
    -resume

# ----------------------------------------------------------
#  5. Surface the report
# ----------------------------------------------------------
REPORT="results_pbmc1k/pbmc1k/report/pbmc1k.report.html"
if [[ -f "$REPORT" ]]; then
    echo
    echo "============================================================"
    echo "  SUCCESS"
    echo "  Per-sample report : $REPORT"
    echo "  MultiQC aggregate : results_pbmc1k/multiqc/multiqc_report.html"
    echo "============================================================"
else
    echo "Report missing — check Nextflow logs above." >&2
    exit 1
fi

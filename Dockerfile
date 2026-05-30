# nf-scrnaseq runtime container
#
# Single image that bundles every binary referenced by the pipeline:
#   - FastQC          (Stage 1)
#   - kb-python       (Stage 2: kallisto|bustools)
#   - scanpy, scrublet, anndata, leidenalg, umap-learn   (Stages 3-5)
#   - MultiQC         (Stage 6)
#
# Build:
#   docker build -t nf-scrnaseq:0.1.0 .
#
# The image targets linux/amd64 for cluster portability. Pin a
# concrete micromamba base for reproducible builds.

FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="nf-scrnaseq"
LABEL org.opencontainers.image.description="scRNA-seq pipeline runtime"
LABEL org.opencontainers.image.licenses="GPL-3.0-or-later"
LABEL org.opencontainers.image.source="https://github.com/example/nf-scrnaseq"

USER root

# Procps gives us `ps`, which Nextflow uses to harvest per-task metrics.
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml

RUN micromamba install -y -n base -f /tmp/env.yml \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
ENV MPLBACKEND=Agg
ENV PYTHONDONTWRITEBYTECODE=1

# Sanity check: every binary the pipeline calls must exist on $PATH
RUN fastqc --version \
    && kb --version \
    && multiqc --version \
    && python -c "import scanpy, scrublet, anndata, leidenalg, umap; print('ok')"

WORKDIR /work

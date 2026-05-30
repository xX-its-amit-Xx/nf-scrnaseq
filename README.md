# nf-scrnaseq

**A reproducible, restartable Nextflow DSL2 pipeline for droplet-based
single-cell RNA-seq.** Takes raw FASTQs → QC'd, clustered AnnData object +
HTML report per sample, plus an aggregated MultiQC report across the cohort.

```
FASTQs  ──▶  FastQC          ─┐
                              ├──▶  MultiQC report
        ──▶  kb-python count ─┤
                  │           ─┘
                  ▼
              Scanpy QC + Scrublet
                  │
                  ▼
        Normalize / HVG / PCA / Leiden / UMAP
                  │
                  ▼
            HTML report per sample
```

---

## What it is

A modular DSL2 pipeline that runs every sample independently and in parallel,
caches every step (`-resume` works), and publishes a clean per-sample
directory layout. Containerized end-to-end via a single Docker image that
bundles FastQC, kb-python, scanpy + scrublet, and MultiQC.

## Who it's for

- Computational biologists processing many samples × tissues uniformly
- Core facilities running a standard scRNA-seq pre-processing service
- Method developers who need a reliable baseline to compare against

## Tech stack & rationale

| Stage              | Tool                | Why                                                              |
| ------------------ | ------------------- | ---------------------------------------------------------------- |
| QC                 | FastQC              | De facto standard, MultiQC-native                                |
| Align + count      | **kb-python**       | Pseudoalignment, ~50 MB index, fast; STARsolo drop-in noted below |
| Per-cell QC        | scanpy              | Canonical Python scRNA-seq stack                                 |
| Doublet detection  | scrublet            | Lightweight, no GPU                                              |
| Cluster + embed    | scanpy (Leiden+UMAP) | Reproducible, well-documented                                    |
| Aggregate report   | MultiQC             | One HTML for the whole cohort                                    |

**Why kb-python instead of STARsolo?** Footprint. The kallisto index for a
human transcriptome is ~3 GB vs STAR's ~30 GB. For routine droplet scRNA-seq
QC + clustering this is a wash on downstream results. If you need
allele-aware quantification, intron counting, or velocyto-style spliced/unspliced
matrices, swap `modules/align_count.nf` for a STARsolo process — the channels
into and out of it are aligner-agnostic.

## Quickstart

```bash
# 1. Build the container (once)
docker build -t nf-scrnaseq:0.1.0 .

# 2. Prepare a samplesheet (see "Samplesheet format" below)
# 3. Run it
nextflow run . \
    --samplesheet samples.csv \
    --kb_index /refs/human/index.idx \
    --t2g     /refs/human/t2g.txt \
    -profile docker
```

For a smoke test on synthetic data:

```bash
python test/prepare_test_data.py           # one-time, ~30 s
nextflow run . -profile test,docker        # ~3 min end-to-end
```

## Samplesheet format

A CSV with one row per sample:

```csv
sample_id,fastq_1,fastq_2,tissue
pbmc_donor1,/data/pbmc_d1_R1.fastq.gz,/data/pbmc_d1_R2.fastq.gz,pbmc
pbmc_donor2,/data/pbmc_d2_R1.fastq.gz,/data/pbmc_d2_R2.fastq.gz,pbmc
lung_donor1,/data/lung_d1_R1.fastq.gz,/data/lung_d1_R2.fastq.gz,lung
```

- `sample_id` — unique per row, used as the publish directory key
- `fastq_1` / `fastq_2` — absolute paths to gzipped FASTQs (10x v2/v3 chemistry)
- `tissue` — free-text label, propagated into the report

Multiple lanes per sample: concatenate beforehand (`zcat L00*_R1*.fastq.gz | gzip > merged_R1.fastq.gz`)
or extend the samplesheet schema and the parser in `main.nf`.

## Profiles

| Profile     | Executor | Container   | Use case                          |
| ----------- | -------- | ----------- | --------------------------------- |
| `standard`  | local    | none        | All tools installed natively      |
| `docker`    | local    | Docker      | Local laptop / dev workstation    |
| `slurm`     | SLURM    | Singularity | HPC cluster                       |
| `test`      | local    | (any)       | Synthetic smoke test, ~3 min      |

Profiles compose: `-profile test,docker` runs the test data in the container.

## Parameter reference

| Parameter         | Default       | Meaning                                                    |
| ----------------- | ------------- | ---------------------------------------------------------- |
| `--samplesheet`   | *required*    | CSV described above                                        |
| `--kb_index`      | *required*    | Pre-built kallisto index (`.idx`)                          |
| `--t2g`           | *required*    | Transcript → gene mapping from `kb ref`                    |
| `--chemistry`     | `10xv3`       | Passed to `kb count -x` (e.g. `10xv2`, `DROPSEQ`)          |
| `--min_genes`     | `200`         | Drop cells with fewer detected genes                       |
| `--min_cells`     | `3`           | Drop genes detected in fewer than this many cells          |
| `--max_mito_pct`  | `20.0`        | Drop cells above this mitochondrial UMI percent            |
| `--doublet_thresh`| `0.25`        | Scrublet score above which a cell is filtered              |
| `--n_top_genes`   | `2000`        | Highly-variable gene count                                 |
| `--n_pcs`         | `50`          | PCA components for the neighborhood graph                  |
| `--n_neighbors`   | `15`          | k for kNN graph                                            |
| `--resolution`    | `0.8`         | Leiden resolution                                          |
| `--outdir`        | `results`     | Where to publish per-sample outputs                        |
| `--publish_mode`  | `copy`        | `copy`, `move`, `link`, `symlink`                          |

All params are overridable on the CLI (`--min_genes 500`) or via a
`-params-file params.yaml`.

## Output layout

```
results/
├── pbmc_donor1/
│   ├── fastqc/                  pbmc_donor1_R1_fastqc.html ...
│   ├── align_count/             counts_unfiltered/adata.h5ad, run_info.json
│   ├── scanpy_qc/               pbmc_donor1.qc.h5ad, pbmc_donor1.qc.json, *_plots/
│   ├── normalize_cluster/       pbmc_donor1.clustered.h5ad, *_plots/
│   └── report/                  pbmc_donor1.report.html
├── pbmc_donor2/ ...
├── multiqc/                     multiqc_report.html, multiqc_data/
└── pipeline_info/               timeline.html, report.html, trace.txt, dag.svg
```

## Resuming and re-running

Every process is cache-keyed on its inputs and parameters, so

```bash
nextflow run . --samplesheet samples.csv -resume
```

skips anything already computed. Tuning `--resolution` re-runs only
`NORMALIZE_CLUSTER` and `REPORT`; tightening `--max_mito_pct` re-runs
`SCANPY_QC` onward; swapping reference triggers everything from
`ALIGN_COUNT` down.

## Building a reference for kb-python

```bash
kb ref \
    -i  index.idx \
    -g  t2g.txt \
    -f1 cdna.fa \
    Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
    Homo_sapiens.GRCh38.110.gtf.gz
```

Cache the result (`index.idx`, `t2g.txt`) somewhere shared — every pipeline
run re-uses it.

## Cookbook

See [COOKBOOK.md](COOKBOOK.md) for end-to-end real-world recipes (10x v3
PBMC, multi-tissue cohort, HPC submission, parameter sweeps).

## License

GNU GPL v3.0 — see [LICENSE](LICENSE).

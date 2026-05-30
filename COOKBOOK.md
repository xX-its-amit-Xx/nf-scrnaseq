# nf-scrnaseq cookbook

Real-world recipes that run end-to-end. Each recipe lists the data, the
exact commands, the expected runtime, and what to look for in the output.

---

## Recipe 1 — 30-second smoke test on synthetic data

**When to use:** verifying your install, sanity-checking after a code change,
CI.

```bash
git clone https://github.com/example/nf-scrnaseq && cd nf-scrnaseq
docker build -t nf-scrnaseq:0.1.0 .
python test/prepare_test_data.py            # build the synthetic dataset
nextflow run . -profile test,docker
```

**Expected:** ~3 min total. `results_test/sampleA/report/sampleA.report.html`
exists, every Leiden cluster has at least one cell, `multiqc_report.html`
shows FastQC sections for 4 FASTQs.

**Pass/fail criteria:**

```bash
# 1. Every per-sample report exists
test -f results_test/sampleA/report/sampleA.report.html
test -f results_test/sampleB/report/sampleB.report.html

# 2. MultiQC saw both samples × both reads
grep -c "fastqc" results_test/multiqc/multiqc_data/multiqc_sources.txt
# -> 4

# 3. Clustered h5ad has a 'leiden' column
python -c "import scanpy as sc; a = sc.read_h5ad('results_test/sampleA/normalize_cluster/sampleA.clustered.h5ad'); assert 'leiden' in a.obs.columns; print(a.obs['leiden'].nunique(), 'clusters')"
```

---

## Recipe 2 — One 10x v3 PBMC sample, local laptop

**Data:** 10x Genomics public `pbmc_1k_v3` (~1.5 GB FASTQs, available from
[support.10xgenomics.com/single-cell-gene-expression/datasets](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3)).

```bash
mkdir -p data/pbmc1k refs/human
cd data/pbmc1k
curl -O https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
tar xf pbmc_1k_v3_fastqs.tar

# Build the reference once (~25 min, 32 GB RAM)
cd ../../refs/human
curl -O http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
curl -O http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
docker run --rm -v "$PWD":/work -w /work nf-scrnaseq:0.1.0 \
    kb ref -i index.idx -g t2g.txt -f1 cdna.fa \
           Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
           Homo_sapiens.GRCh38.110.gtf.gz

# Samplesheet
cd ../..
cat > samples.csv <<EOF
sample_id,fastq_1,fastq_2,tissue
pbmc1k,$(pwd)/data/pbmc1k/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,$(pwd)/data/pbmc1k/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz,pbmc
EOF

nextflow run . \
    --samplesheet  samples.csv \
    --kb_index     $(pwd)/refs/human/index.idx \
    --t2g          $(pwd)/refs/human/t2g.txt \
    --chemistry    10xv3 \
    -profile       docker
```

**Expected runtime:** ~45 min on a 16-core / 32 GB workstation
(`ALIGN_COUNT` is the long pole at ~30 min).

**What to check:**

- `results/pbmc1k/scanpy_qc/pbmc1k.qc.json` → `n_cells_post` between 900
  and 1200 (10x advertises ~1000)
- `results/pbmc1k/report/pbmc1k.report.html` → UMAP shows distinct lymphoid
  vs myeloid lobes; expect 8–12 Leiden clusters at default resolution
- `results/multiqc/multiqc_report.html` → R1 per-base quality drops near
  3' (this is the expected 10x barcode read profile, not a problem)

---

## Recipe 3 — Multi-tissue cohort with merged lanes

**Data:** four samples across two tissues, two donors each. Two of them
have been sequenced across two lanes.

```bash
# Pre-merge the lanes (Nextflow could do this in a process; for clarity
# we do it once up front).
for s in pbmc_d1 lung_d1; do
    zcat data/${s}_L00{1,2}_R1_001.fastq.gz | gzip > data/merged/${s}_R1.fastq.gz
    zcat data/${s}_L00{1,2}_R2_001.fastq.gz | gzip > data/merged/${s}_R2.fastq.gz
done
# Single-lane samples just get symlinked
ln -sf $(pwd)/data/pbmc_d2_L001_R1_001.fastq.gz data/merged/pbmc_d2_R1.fastq.gz
ln -sf $(pwd)/data/pbmc_d2_L001_R2_001.fastq.gz data/merged/pbmc_d2_R2.fastq.gz
# ... same for lung_d2

cat > cohort.csv <<EOF
sample_id,fastq_1,fastq_2,tissue
pbmc_d1,$(pwd)/data/merged/pbmc_d1_R1.fastq.gz,$(pwd)/data/merged/pbmc_d1_R2.fastq.gz,pbmc
pbmc_d2,$(pwd)/data/merged/pbmc_d2_R1.fastq.gz,$(pwd)/data/merged/pbmc_d2_R2.fastq.gz,pbmc
lung_d1,$(pwd)/data/merged/lung_d1_R1.fastq.gz,$(pwd)/data/merged/lung_d1_R2.fastq.gz,lung
lung_d2,$(pwd)/data/merged/lung_d2_R1.fastq.gz,$(pwd)/data/merged/lung_d2_R2.fastq.gz,lung
EOF

# Lung samples tend to have higher mito% — relax that filter via params file
cat > cohort.params.yaml <<EOF
samplesheet:   cohort.csv
kb_index:      /refs/human/index.idx
t2g:           /refs/human/t2g.txt
max_mito_pct:  25.0
min_genes:     300
resolution:    1.0
outdir:        cohort_results
EOF

nextflow run . -params-file cohort.params.yaml -profile docker
```

**Expected:** all four samples are processed in parallel up to the
`ALIGN_COUNT` slot limit. Per-sample reports land in
`cohort_results/<sample_id>/`; the MultiQC report shows all four FastQC
groups side-by-side.

**Tip:** because tissue is a column in the samplesheet, you can recover
it from the AnnData later for cross-tissue integration:

```python
import anndata as ad
adatas = {sid: ad.read_h5ad(f'cohort_results/{sid}/normalize_cluster/{sid}.clustered.h5ad')
          for sid in ['pbmc_d1','pbmc_d2','lung_d1','lung_d2']}
combined = ad.concat(adatas, label='sample_id', join='outer')
# combined.obs['tissue'] propagates from each AnnData's per-cell annotation
```

---

## Recipe 4 — SLURM cluster, 24-sample cohort

```bash
# Stage data on shared scratch
sbatch --wrap "rsync -av /archive/seq/run42/ /scratch/$USER/run42/"

# samplesheet generated from the manifest
awk -F, 'NR>1 {
   printf "%s,/scratch/'$USER'/run42/%s_R1.fastq.gz,/scratch/'$USER'/run42/%s_R2.fastq.gz,%s\n",
          $1,$1,$1,$2
}' manifest.csv > run42.csv
sed -i '1i sample_id,fastq_1,fastq_2,tissue' run42.csv

# Launch — Nextflow itself runs on the head node; per-process jobs go to SLURM
nextflow run /opt/pipelines/nf-scrnaseq \
    --samplesheet run42.csv \
    --kb_index    /shared/refs/human/index.idx \
    --t2g         /shared/refs/human/t2g.txt \
    --outdir      /scratch/$USER/run42_results \
    -profile      slurm \
    -with-tower    # if you have Nextflow Tower configured
```

**Capacity notes:** `conf/slurm.config` caps `executor.queueSize = 50` — for
a 24-sample cohort that means all alignments queue together. Bump it if
your fair-share allows; drop it if you're sharing a small partition.

**Re-run a single failed sample:** Nextflow's resume is sample-scoped. If
`lung_d7` failed alignment because the R2 file was truncated, just fix the
file and re-run — the other 23 samples' work is cached and won't recompute.

```bash
nextflow run /opt/pipelines/nf-scrnaseq -params-file run42.yaml -resume
```

---

## Recipe 5 — Parameter sweep for QC thresholds

You're not sure whether `--max_mito_pct 15` or `25` is right for a new
tissue. Run both, compare, pick.

```bash
for mito in 10 15 20 25; do
    nextflow run . \
        --samplesheet  samples.csv \
        --kb_index     /refs/human/index.idx \
        --t2g          /refs/human/t2g.txt \
        --max_mito_pct $mito \
        --outdir       results_mito${mito} \
        -profile       docker \
        -resume
done
```

`ALIGN_COUNT` runs once across the four invocations (cache hits) because
the aligner inputs don't change with `--max_mito_pct`. Only
`SCANPY_QC → NORMALIZE_CLUSTER → REPORT` re-runs per sweep value.

Compare:

```bash
for mito in 10 15 20 25; do
    n=$(python -c "import json; print(json.load(open('results_mito${mito}/sampleA/scanpy_qc/sampleA.qc.json'))['n_cells_post'])")
    echo "mito=${mito}%: ${n} cells"
done
```

---

## Recipe 6 — Adding a per-sample tweak via samplesheet

Suppose `lung_d2` was sequenced on v2 chemistry while everything else is v3.
The cleanest extension is to add a `chemistry` column to the samplesheet and
read it in `main.nf`:

```diff
 ch_samples = Channel
     .fromPath(params.samplesheet, checkIfExists: true)
     .splitCsv(header: true)
     .map { row ->
         def meta = [
             id:        row.sample_id,
             tissue:    row.tissue ?: 'unspecified',
+            chemistry: row.chemistry ?: params.chemistry
         ]
         def r1 = file(row.fastq_1, checkIfExists: true)
         def r2 = file(row.fastq_2, checkIfExists: true)
         tuple(meta, r1, r2)
     }
```

then in `modules/align_count.nf`:

```diff
-    -x ${params.chemistry} \\
+    -x ${meta.chemistry} \\
```

Now each sample carries its own chemistry and the pipeline does the right
thing without per-sample command-line gymnastics.

---

## Troubleshooting

| Symptom                                                         | Likely cause                                       | Fix                                                                  |
| --------------------------------------------------------------- | -------------------------------------------------- | -------------------------------------------------------------------- |
| `kb count` fails with `Error: file does not exist`              | Wrong `--kb_index` path                            | Verify with `kb inspect index.idx`                                   |
| `SCANPY_QC` survives 0 cells                                    | Filters too strict for the data quality            | Loosen `--min_genes`, `--max_mito_pct`; re-run with `-resume`        |
| `NORMALIZE_CLUSTER` errors `n_neighbors > n_samples`            | QC dropped too many cells                          | Reduce `--n_neighbors` or fix QC first                               |
| Docker process hangs on Mac                                     | File-sharing layer being slow on bind mounts       | Use `:delegated` mounts or move work dir under a docker-volume       |
| `-resume` re-runs everything anyway                             | A param changed silently (e.g. `outdir` typo)      | Diff `nextflow.config` and `params-file` between runs                |
| FastQC reports "Per base sequence content" failures on R1       | Expected — R1 carries cell barcodes, not biology   | Ignore; check R2 quality instead                                     |

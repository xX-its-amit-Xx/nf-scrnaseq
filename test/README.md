# Test profile

The `-profile test` runs the pipeline against a deterministic synthetic
dataset built by `prepare_test_data.py`. This keeps the git repo small (no
binary FASTQs committed) while still letting reviewers exercise every stage
end-to-end.

## One-time setup

The test data + kallisto index must exist before the first `-profile test`
run. Either run the prep script directly:

```bash
python test/prepare_test_data.py
```

or, if you don't have `kallisto` installed on the host, run it inside the
pipeline's container:

```bash
docker run --rm -v "$PWD":/work -w /work nf-scrnaseq:0.1.0 \
    python test/prepare_test_data.py
```

This produces:

```
test/
├── ref/
│   ├── cdna.fa        # 12 synthetic transcripts, 2 kb each
│   ├── t2g.txt        # transcript -> gene mapping
│   └── index.idx      # kallisto index
├── fastqs/
│   ├── sampleA_R1.fastq.gz   # 16 bp barcode + 12 bp UMI
│   ├── sampleA_R2.fastq.gz   # 90 bp cDNA reads
│   ├── sampleB_R1.fastq.gz
│   └── sampleB_R2.fastq.gz
└── samplesheet.csv
```

## Run the smoke test

```bash
nextflow run . -profile test,docker
```

End-to-end runtime on a laptop is ~3 minutes. The `results_test/` directory
will contain per-sample QC + clustering HTML reports and an aggregated
MultiQC report.

## Why synthetic?

| Option              | Pros                          | Cons                                    |
| ------------------- | ----------------------------- | --------------------------------------- |
| Bundle real FASTQs  | Realistic biology             | 50-500 MB in git, license issues        |
| Download on first run | Realistic biology           | Flaky, network-dependent                |
| **Synthetic (ours)** | Tiny, deterministic, offline | No real biology — clustering is noise |

For this pipeline the test is a *plumbing* check, not a biology check —
every process must run, every channel must connect, every output must
materialize. Synthetic data answers that question reliably.

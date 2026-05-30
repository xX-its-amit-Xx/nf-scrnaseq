# Convenience targets for nf-scrnaseq.
# Designed to work inside the .devcontainer (or any Linux env with
# Docker + Nextflow installed).

IMAGE       ?= nf-scrnaseq:0.1.0
TEST_OUTDIR ?= results_test

.PHONY: image data test pbmc1k clean help

help:
	@echo "Targets:"
	@echo "  make image    Build the runtime container ($(IMAGE))"
	@echo "  make data     Generate synthetic test FASTQs + kallisto index"
	@echo "  make test     End-to-end smoke test against the synthetic dataset"
	@echo "  make pbmc1k   Real-world recipe: 10x pbmc_1k_v3 end-to-end"
	@echo "  make clean    Remove Nextflow work, results, and generated test data"

image:
	docker build -t $(IMAGE) .

data: image
	docker run --rm -v "$$PWD:/work" -w /work $(IMAGE) \
	    python test/prepare_test_data.py

test: data
	nextflow run . -profile test,docker --outdir $(TEST_OUTDIR)
	@echo
	@echo "Reports:"
	@ls -1 $(TEST_OUTDIR)/*/report/*.html 2>/dev/null || true
	@ls -1 $(TEST_OUTDIR)/multiqc/multiqc_report.html 2>/dev/null || true

pbmc1k: image
	bash scripts/run_pbmc1k.sh

clean:
	rm -rf .nextflow* work/ results $(TEST_OUTDIR) \
	       test/ref test/fastqs

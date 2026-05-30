#!/usr/bin/env bash
# Codespaces post-create hook for nf-scrnaseq.
#
# This is intentionally lightweight: only the steps that *must* succeed for
# the Codespace to be usable (Nextflow install) run in the hook itself. The
# heavier optional steps (building the runtime container, generating test
# data) are run best-effort so that even if they fail, the Codespace still
# comes up — you can then debug interactively instead of being dropped into
# a recovery container.
set -uo pipefail

cd "$(dirname "$0")/.."

# ----------------------------------------------------------
# Required: install Nextflow
# ----------------------------------------------------------
echo "::group::Install Nextflow"
if ! command -v nextflow >/dev/null 2>&1; then
    curl -fsSL https://get.nextflow.io | bash
    sudo mv nextflow /usr/local/bin/nextflow
    sudo chmod +x /usr/local/bin/nextflow
fi
nextflow -version || { echo "Nextflow install failed" >&2; exit 1; }
echo "::endgroup::"

# ----------------------------------------------------------
# Best-effort: build the runtime container.
#   `|| true` so a build failure (e.g. transient mirror outage) doesn't
#   block the Codespace from coming up. Re-run with `make image` once
#   you're inside.
# ----------------------------------------------------------
echo "::group::Build nf-scrnaseq container (best-effort)"
docker build -t nf-scrnaseq:0.1.0 . || \
    echo "WARN: docker build failed — re-run with 'make image' once inside the Codespace"
echo "::endgroup::"

# ----------------------------------------------------------
# Best-effort: generate the synthetic test dataset.
# ----------------------------------------------------------
echo "::group::Generate synthetic test dataset (best-effort)"
if docker image inspect nf-scrnaseq:0.1.0 >/dev/null 2>&1; then
    docker run --rm -v "$PWD:/work" -w /work nf-scrnaseq:0.1.0 \
        python test/prepare_test_data.py \
        || echo "WARN: test data prep failed — re-run with 'make data'"
else
    echo "Skipping: nf-scrnaseq:0.1.0 image not built"
fi
echo "::endgroup::"

cat <<'EOF'

============================================================
  nf-scrnaseq Codespace is ready.

  Smoke test (~3 min):
      make test                  # or: nextflow run . -profile test,docker

  Real-world recipe (10x PBMC 1k, ~30 min):
      make pbmc1k                # or: bash scripts/run_pbmc1k.sh

  See README.md and COOKBOOK.md.
============================================================
EOF

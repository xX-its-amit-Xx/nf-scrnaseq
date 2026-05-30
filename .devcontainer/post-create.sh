#!/usr/bin/env bash
# Codespaces post-create hook for nf-scrnaseq.
#
# Runs once when the dev container is first built. Idempotent — safe to
# re-run if you blow away the install.
set -euo pipefail

cd "$(dirname "$0")/.."

echo "::group::Install Nextflow"
if ! command -v nextflow >/dev/null 2>&1; then
    curl -fsSL https://get.nextflow.io | bash
    sudo mv nextflow /usr/local/bin/nextflow
    sudo chmod +x /usr/local/bin/nextflow
fi
nextflow -version
echo "::endgroup::"

echo "::group::Build nf-scrnaseq container"
docker build -t nf-scrnaseq:0.1.0 .
echo "::endgroup::"

echo "::group::Generate synthetic test dataset"
docker run --rm -v "$PWD:/work" -w /work nf-scrnaseq:0.1.0 \
    python test/prepare_test_data.py
echo "::endgroup::"

cat <<'EOF'

============================================================
  nf-scrnaseq Codespace is ready.

  Smoke test (end-to-end, ~3 min):
      nextflow run . -profile test,docker

  Real run example:
      nextflow run . \
          --samplesheet samples.csv \
          --kb_index    /refs/index.idx \
          --t2g         /refs/t2g.txt \
          -profile      docker

  See README.md and COOKBOOK.md.
============================================================
EOF

#!/usr/bin/env python3
"""
prepare_test_data.py — generate a deterministic synthetic dataset for the
nf-scrnaseq `-profile test` smoke test.

What it produces (all under test/):
    ref/cdna.fa                 — 12 synthetic "transcripts", 2 kb each
    ref/index.idx               — kallisto index built from cdna.fa
    ref/t2g.txt                 — transcript -> gene mapping
    fastqs/sampleA_R1.fastq.gz  — 10x v3 reads: 16 bp barcode + 12 bp UMI
    fastqs/sampleA_R2.fastq.gz  — cDNA reads (90 bp)
    fastqs/sampleB_R1.fastq.gz
    fastqs/sampleB_R2.fastq.gz
    samplesheet.csv             — points at the four files above

Requirements:
    - kallisto on $PATH (provided by the Docker image; on the host you can
      use `conda install -c bioconda kallisto`)

Determinism: all randomness is seeded so two runs produce identical FASTQs.
"""
from __future__ import annotations

import argparse
import gzip
import random
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable

# ---------------------------------------------------------------------------
#  Reference generation
# ---------------------------------------------------------------------------

N_TRANSCRIPTS = 12
TRANSCRIPT_LEN = 2000
ALPHABET = "ACGT"


def synthesize_transcripts(rng: random.Random) -> dict[str, str]:
    """Return {transcript_id: sequence} with each transcript a random string.

    Each transcript begins with a unique 'tag' prefix so pseudoalignment of
    the simulated reads (which include the tag) is unambiguous.
    """
    transcripts = {}
    for i in range(N_TRANSCRIPTS):
        tid = f"ENST_TEST_{i:03d}"
        tag = f"TTGCA{i:03d}AAGCAT".replace(" ", "")  # transcript-specific prefix
        body = "".join(rng.choices(ALPHABET, k=TRANSCRIPT_LEN - len(tag)))
        transcripts[tid] = tag + body
    return transcripts


def write_fasta(path: Path, seqs: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for name, seq in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def write_t2g(path: Path, transcripts: dict[str, str]) -> None:
    """t2g format expected by kb-python: transcript<TAB>gene<TAB>gene_name."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for tid in transcripts:
            gid = tid.replace("ENST", "ENSG")
            fh.write(f"{tid}\t{gid}\t{gid}\n")


# ---------------------------------------------------------------------------
#  FASTQ simulation (10x v3 chemistry)
# ---------------------------------------------------------------------------

CB_LEN = 16
UMI_LEN = 12
R2_LEN = 90


def random_barcode(rng: random.Random, length: int) -> str:
    return "".join(rng.choices(ALPHABET, k=length))


def load_10x_v3_whitelist() -> list[str] | None:
    """Return a list of valid 10x v3 cell barcodes, or None if unavailable.

    The whitelist ships with ngs_tools (a kb-python dep). Random 16-mers
    would mostly be rejected by kb's barcode-correction step (the v3
    whitelist covers ~6.8M of 4^16 possible 16-mers, so the hit rate is
    ~0.16%). Drawing barcodes from the whitelist makes the synthetic
    dataset survive kb count and produce a non-degenerate test report.
    """
    try:
        import ngs_tools  # noqa: F401  (we just need its install dir)
    except ImportError:
        return None
    ngs_dir = Path(__import__("ngs_tools").__file__).parent
    candidates = list(ngs_dir.glob("**/10x_version3_whitelist.txt*"))
    if not candidates:
        return None
    path = candidates[0]
    opener: type = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [line.strip() for line in fh if line.strip()]


def sample_barcodes(rng: random.Random, n: int) -> list[str]:
    """Prefer real 10x v3 whitelist barcodes; fall back to random 16-mers."""
    whitelist = load_10x_v3_whitelist()
    if whitelist is None:
        print(
            "[prepare_test_data] kb_python whitelist not found; using random barcodes",
            file=sys.stderr,
        )
        return [random_barcode(rng, CB_LEN) for _ in range(n)]
    return rng.sample(whitelist, k=n)


def fastq_record(read_id: str, seq: str) -> str:
    qual = "I" * len(seq)  # Phred 40, plenty for kallisto pseudoalignment
    return f"@{read_id}\n{seq}\n+\n{qual}\n"


def simulate_sample(
    sample_id: str,
    transcripts: dict[str, str],
    n_cells: int,
    reads_per_cell: int,
    out_r1: Path,
    out_r2: Path,
    rng: random.Random,
) -> None:
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    transcript_ids = list(transcripts.keys())
    barcodes = sample_barcodes(rng, n_cells)

    with gzip.open(out_r1, "wt") as r1_fh, gzip.open(out_r2, "wt") as r2_fh:
        for cell_idx, barcode in enumerate(barcodes):
            # Each "cell" expresses a random subset of transcripts so
            # post-clustering UMAP shows real structure on this toy dataset.
            # Keep the expressed set wider than the test profile's min_genes
            # so cells survive QC end-to-end.
            expressed = rng.sample(
                transcript_ids,
                k=rng.randint(6, min(10, len(transcript_ids)))
            )
            for read_idx in range(reads_per_cell):
                tid = rng.choice(expressed)
                seq = transcripts[tid]
                start = rng.randint(0, len(seq) - R2_LEN)
                r2 = seq[start:start + R2_LEN]
                umi = random_barcode(rng, UMI_LEN)
                rid = f"{sample_id}:cell{cell_idx}:read{read_idx}"
                r1_fh.write(fastq_record(rid, barcode + umi))
                r2_fh.write(fastq_record(rid, r2))


# ---------------------------------------------------------------------------
#  kallisto index
# ---------------------------------------------------------------------------

def build_kallisto_index(cdna_fa: Path, out_idx: Path) -> None:
    if shutil.which("kallisto") is None:
        sys.exit(
            "kallisto not found on $PATH. Install it (e.g. "
            "`conda install -c bioconda kallisto`) or run this inside the "
            "nf-scrnaseq Docker image."
        )
    out_idx.parent.mkdir(parents=True, exist_ok=True)
    # k=31 is the kallisto default; transcripts are 2kb so this is fine.
    subprocess.run(
        ["kallisto", "index", "-i", str(out_idx), str(cdna_fa)],
        check=True,
    )


# ---------------------------------------------------------------------------
#  Samplesheet
# ---------------------------------------------------------------------------

def write_samplesheet(path: Path, samples: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("sample_id,fastq_1,fastq_2,tissue\n")
        for s in samples:
            fh.write(
                f"{s['id']},{s['r1']},{s['r2']},{s['tissue']}\n"
            )


# ---------------------------------------------------------------------------
#  Main
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", type=Path, default=Path(__file__).parent)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--n-cells", type=int, default=200)
    p.add_argument("--reads-per-cell", type=int, default=80)
    p.add_argument(
        "--skip-index",
        action="store_true",
        help="Skip kallisto index build (useful when iterating on the simulator)",
    )
    args = p.parse_args()

    rng = random.Random(args.seed)

    print(f"[prepare_test_data] writing into {args.out_dir.resolve()}")

    transcripts = synthesize_transcripts(rng)
    cdna_fa = args.out_dir / "ref" / "cdna.fa"
    t2g     = args.out_dir / "ref" / "t2g.txt"
    idx     = args.out_dir / "ref" / "index.idx"
    write_fasta(cdna_fa, transcripts)
    write_t2g(t2g, transcripts)
    print(f"  ref:  {cdna_fa.name}, {t2g.name}")

    if not args.skip_index:
        build_kallisto_index(cdna_fa, idx)
        print(f"  ref:  {idx.name}")

    # Write paths relative to the project root so the samplesheet works
    # whether you invoke the script from the host or from inside the
    # container (where the resolved absolute path would be /work/...).
    project_root = args.out_dir.parent

    samples = []
    for sid, tissue in [("sampleA", "pbmc"), ("sampleB", "lung")]:
        r1 = args.out_dir / "fastqs" / f"{sid}_R1.fastq.gz"
        r2 = args.out_dir / "fastqs" / f"{sid}_R2.fastq.gz"
        simulate_sample(
            sid, transcripts, args.n_cells, args.reads_per_cell,
            r1, r2, rng,
        )
        samples.append({
            "id": sid,
            "tissue": tissue,
            "r1": r1.resolve().relative_to(project_root.resolve()).as_posix(),
            "r2": r2.resolve().relative_to(project_root.resolve()).as_posix(),
        })
        print(f"  fastq: {r1.name}, {r2.name}")

    samplesheet = args.out_dir / "samplesheet.csv"
    write_samplesheet(samplesheet, samples)
    print(f"  samplesheet: {samplesheet.name}")

    print("\nReady. Run:")
    print("  nextflow run . -profile test,docker")
    return 0


if __name__ == "__main__":
    sys.exit(main())

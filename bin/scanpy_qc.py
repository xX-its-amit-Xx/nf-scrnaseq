#!/usr/bin/env python3
"""
scanpy_qc.py — load a kb-python count matrix, compute per-cell QC metrics,
run Scrublet for doublet detection, and apply configurable filters.

Inputs:
    --counts            kb-python counts_unfiltered/ directory (contains adata.h5ad
                        and/or .mtx + barcodes.txt + genes.txt)
Outputs:
    --out-h5ad          QC'd AnnData (raw counts retained in .layers['counts'])
    --out-json          QC summary metrics (used by REPORT and MultiQC)
    --out-plotdir       Per-sample diagnostic PNGs (violin, scatter, knee)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

sc.settings.verbosity = 1


def load_counts(counts_dir: Path) -> ad.AnnData:
    """Load a kb-python output directory into AnnData.

    kb-python emits adata.h5ad when run with --h5ad. If that file is missing,
    fall back to assembling from the MTX + barcodes + genes triple.
    """
    h5ad = counts_dir / "adata.h5ad"
    if h5ad.exists():
        return sc.read_h5ad(h5ad)

    mtx = counts_dir / "cells_x_genes.mtx"
    if not mtx.exists():
        raise FileNotFoundError(
            f"No adata.h5ad or cells_x_genes.mtx in {counts_dir}. "
            f"Contents: {list(counts_dir.iterdir())}"
        )

    adata = sc.read_mtx(mtx).T if _needs_transpose(mtx) else sc.read_mtx(mtx)
    barcodes = (counts_dir / "cells_x_genes.barcodes.txt").read_text().splitlines()
    genes = (counts_dir / "cells_x_genes.genes.txt").read_text().splitlines()
    adata.obs_names = barcodes
    adata.var_names = genes
    return adata


def _needs_transpose(mtx_path: Path) -> bool:
    # kb-python writes cells_x_genes.mtx with cells as rows; nothing to do.
    return False


def run_scrublet(adata: ad.AnnData, threshold: float) -> ad.AnnData:
    """Compute Scrublet doublet scores; tolerate small-sample failures."""
    try:
        import scrublet as scr
    except ImportError:
        print("[scanpy_qc] scrublet not installed; skipping", file=sys.stderr)
        adata.obs["doublet_score"] = 0.0
        adata.obs["predicted_doublet"] = False
        return adata

    counts = adata.layers.get("counts", adata.X)
    try:
        scrub = scr.Scrublet(counts)
        scores, predicted = scrub.scrub_doublets(verbose=False)
        adata.obs["doublet_score"] = scores
        adata.obs["predicted_doublet"] = (
            predicted if predicted is not None else scores > threshold
        )
    except Exception as e:  # noqa: BLE001 — Scrublet can crash on tiny matrices
        print(f"[scanpy_qc] Scrublet failed ({e}); zero-filling", file=sys.stderr)
        adata.obs["doublet_score"] = 0.0
        adata.obs["predicted_doublet"] = False
    return adata


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--counts", required=True, type=Path)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--tissue", default="unspecified")
    p.add_argument("--min-genes", type=int, default=200)
    p.add_argument("--min-cells", type=int, default=3)
    p.add_argument("--max-mito-pct", type=float, default=20.0)
    p.add_argument("--doublet-thresh", type=float, default=0.25)
    p.add_argument("--out-h5ad", required=True, type=Path)
    p.add_argument("--out-json", required=True, type=Path)
    p.add_argument("--out-plotdir", required=True, type=Path)
    args = p.parse_args()

    args.out_plotdir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = args.out_plotdir

    adata = load_counts(args.counts)
    adata.var_names_make_unique()
    adata.obs["sample_id"] = args.sample_id
    adata.obs["tissue"] = args.tissue
    adata.layers["counts"] = adata.X.copy()

    n_cells_pre = adata.n_obs
    n_genes_pre = adata.n_vars

    # Mitochondrial gene flag — works for human (MT-) and mouse (mt-) symbols
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # Pre-filter diagnostic plots — tolerated to fail on tiny / degenerate
    # data (seaborn 0.13+ catplot can throw IndexError when a group has 0
    # rows). The QC json and h5ad are the load-bearing outputs.
    try:
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show=False,
            save=f"_{args.sample_id}_prefilter.png",
        )
    except Exception as e:  # noqa: BLE001
        print(f"[scanpy_qc] pre-filter violin failed ({e}); skipping", file=sys.stderr)

    adata = run_scrublet(adata, args.doublet_thresh)

    # Apply filters
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    adata = adata[adata.obs["pct_counts_mt"] <= args.max_mito_pct].copy()
    adata = adata[adata.obs["doublet_score"] <= args.doublet_thresh].copy()

    # Post-filter diagnostic plot — same tolerance as the pre-filter one.
    if adata.n_obs > 0:
        try:
            sc.pl.violin(
                adata,
                ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                jitter=0.4,
                multi_panel=True,
                show=False,
                save=f"_{args.sample_id}_postfilter.png",
            )
        except Exception as e:  # noqa: BLE001
            print(f"[scanpy_qc] post-filter violin failed ({e}); skipping", file=sys.stderr)

    metrics = {
        "sample_id": args.sample_id,
        "tissue": args.tissue,
        "n_cells_pre": int(n_cells_pre),
        "n_genes_pre": int(n_genes_pre),
        "n_cells_post": int(adata.n_obs),
        "n_genes_post": int(adata.n_vars),
        "median_genes_per_cell": float(np.median(adata.obs["n_genes_by_counts"]))
        if adata.n_obs
        else 0.0,
        "median_counts_per_cell": float(np.median(adata.obs["total_counts"]))
        if adata.n_obs
        else 0.0,
        "median_pct_mito": float(np.median(adata.obs["pct_counts_mt"]))
        if adata.n_obs
        else 0.0,
        "n_predicted_doublets": int(adata.obs.get("predicted_doublet", []).sum())
        if adata.n_obs
        else 0,
        "filters": {
            "min_genes": args.min_genes,
            "min_cells": args.min_cells,
            "max_mito_pct": args.max_mito_pct,
            "doublet_thresh": args.doublet_thresh,
        },
    }

    adata.write_h5ad(args.out_h5ad)
    args.out_json.write_text(json.dumps(metrics, indent=2))
    print(f"[scanpy_qc] {args.sample_id}: {n_cells_pre} -> {adata.n_obs} cells")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
normalize_cluster.py — normalize counts, find HVGs, scale, PCA, neighbors,
Leiden clustering, and UMAP embedding.

Workflow follows the Scanpy recommended pipeline:
    counts -> normalize_total(1e4) -> log1p -> HVG -> scale -> PCA ->
    neighbors -> leiden -> UMAP
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import scanpy as sc

sc.settings.verbosity = 1


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--in-h5ad", required=True, type=Path)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--n-top-genes", type=int, default=2000)
    p.add_argument("--n-pcs", type=int, default=50)
    p.add_argument("--n-neighbors", type=int, default=15)
    p.add_argument("--resolution", type=float, default=0.8)
    p.add_argument("--out-h5ad", required=True, type=Path)
    p.add_argument("--out-plotdir", required=True, type=Path)
    args = p.parse_args()

    args.out_plotdir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = args.out_plotdir

    adata = sc.read_h5ad(args.in_h5ad)
    if adata.n_obs < 10:
        print(
            f"[normalize_cluster] Only {adata.n_obs} cells survived QC; "
            f"writing input through unchanged",
            file=sys.stderr,
        )
        adata.write_h5ad(args.out_h5ad)
        return 0

    # Stash raw counts before in-place normalization
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    # HVG selection — clamp to available gene count for the tiny test profile
    n_top = min(args.n_top_genes, adata.n_vars - 1)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat")
    sc.pl.highly_variable_genes(
        adata, show=False, save=f"_{args.sample_id}_hvg.png"
    )

    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)

    n_pcs = min(args.n_pcs, max(adata.n_obs, adata.n_vars) - 1, 50)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
    sc.pl.pca_variance_ratio(
        adata, n_pcs=n_pcs, show=False, save=f"_{args.sample_id}_pca_var.png"
    )

    n_neighbors = min(args.n_neighbors, max(adata.n_obs - 1, 2))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.leiden(adata, resolution=args.resolution, key_added="leiden")
    sc.tl.umap(adata)

    sc.pl.umap(
        adata,
        color=["leiden"],
        show=False,
        save=f"_{args.sample_id}_umap_leiden.png",
    )
    if "pct_counts_mt" in adata.obs:
        sc.pl.umap(
            adata,
            color=["pct_counts_mt", "n_genes_by_counts"],
            show=False,
            save=f"_{args.sample_id}_umap_qc.png",
        )

    adata.write_h5ad(args.out_h5ad)
    print(
        f"[normalize_cluster] {args.sample_id}: "
        f"{adata.n_obs} cells, {adata.obs['leiden'].nunique()} clusters"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

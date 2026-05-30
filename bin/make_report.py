#!/usr/bin/env python3
"""
make_report.py — render a self-contained HTML QC + clustering report.

Reads:
    - a clustered AnnData (.h5ad)
    - a QC metrics JSON (from scanpy_qc.py)
    - an HTML template with ${placeholder} tokens

Writes a single-file HTML report with embedded PNGs (base64) so the
artifact is portable and survives publishDir copies.
"""
from __future__ import annotations

import argparse
import base64
import io
import json
import sys
from pathlib import Path
from string import Template

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


def fig_to_b64(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=110)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def render_metrics_table(metrics: dict) -> str:
    rows = [
        ("Sample ID", metrics.get("sample_id", "-")),
        ("Tissue", metrics.get("tissue", "-")),
        ("Cells (pre-QC)", f"{metrics.get('n_cells_pre', 0):,}"),
        ("Cells (post-QC)", f"{metrics.get('n_cells_post', 0):,}"),
        ("Genes (post-QC)", f"{metrics.get('n_genes_post', 0):,}"),
        ("Median genes/cell", f"{metrics.get('median_genes_per_cell', 0):.0f}"),
        ("Median counts/cell", f"{metrics.get('median_counts_per_cell', 0):.0f}"),
        ("Median % mito", f"{metrics.get('median_pct_mito', 0):.2f}"),
        ("Predicted doublets", f"{metrics.get('n_predicted_doublets', 0):,}"),
    ]
    return "\n".join(
        f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in rows
    )


def render_filters_table(filters: dict) -> str:
    return "\n".join(
        f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in filters.items()
    )


def umap_panel(adata: ad.AnnData, color: str, title: str) -> str:
    # Skip silently when this AnnData never got embedded — e.g. a sample
    # that fell below the cell-count cutoff in NORMALIZE_CLUSTER.
    if "X_umap" not in adata.obsm or color not in adata.obs.columns:
        return ""
    try:
        fig = sc.pl.umap(
            adata, color=color, return_fig=True, show=False, title=title
        )
    except Exception as e:  # noqa: BLE001
        return (
            f'<div class="panel"><h3>{title}</h3>'
            f'<p><em>plot failed: {e}</em></p></div>'
        )
    return (
        f'<div class="panel"><h3>{title}</h3>'
        f'<img src="data:image/png;base64,{fig_to_b64(fig)}"/></div>'
    )


def qc_violin(adata: ad.AnnData) -> str:
    cols = [
        c for c in ("n_genes_by_counts", "total_counts", "pct_counts_mt")
        if c in adata.obs.columns
    ]
    if not cols:
        return ""
    fig, axes = plt.subplots(1, len(cols), figsize=(4 * len(cols), 4))
    if len(cols) == 1:
        axes = [axes]
    for ax, col in zip(axes, cols):
        # matplotlib.violinplot blows up on empty arrays — fall back to a
        # text label so the rest of the report still renders.
        values = adata.obs[col].dropna().to_numpy()
        if values.size == 0:
            ax.text(0.5, 0.5, "no data", ha="center", va="center")
        else:
            try:
                ax.violinplot(values, showmedians=True)
            except (ValueError, np.linalg.LinAlgError) as e:
                ax.text(0.5, 0.5, f"plot failed:\n{e}", ha="center", va="center", fontsize=8)
        ax.set_title(col)
        ax.set_xticks([])
    fig.tight_layout()
    return (
        f'<div class="panel"><h3>Post-QC distributions</h3>'
        f'<img src="data:image/png;base64,{fig_to_b64(fig)}"/></div>'
    )


def cluster_size_table(adata: ad.AnnData) -> str:
    if "leiden" not in adata.obs.columns:
        return "<p><em>No clustering available.</em></p>"
    counts = adata.obs["leiden"].value_counts().sort_index()
    rows = "\n".join(
        f"<tr><td>{c}</td><td>{n:,}</td></tr>" for c, n in counts.items()
    )
    return (
        "<table><thead><tr><th>Cluster</th><th>Cells</th></tr></thead>"
        f"<tbody>{rows}</tbody></table>"
    )


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True, type=Path)
    p.add_argument("--qc-json", required=True, type=Path)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--tissue", default="unspecified")
    p.add_argument("--template", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    metrics = json.loads(args.qc_json.read_text())
    adata = sc.read_h5ad(args.h5ad)

    panels = []
    panels.append(qc_violin(adata))
    panels.append(umap_panel(adata, "leiden", "UMAP by Leiden cluster"))
    panels.append(umap_panel(adata, "pct_counts_mt", "UMAP by % mito"))
    panels.append(umap_panel(adata, "n_genes_by_counts", "UMAP by gene count"))
    panels.append(umap_panel(adata, "doublet_score", "UMAP by doublet score"))

    template = Template(args.template.read_text())
    html = template.safe_substitute(
        sample_id=args.sample_id,
        tissue=args.tissue,
        n_cells=adata.n_obs,
        n_clusters=int(adata.obs["leiden"].nunique())
        if "leiden" in adata.obs.columns
        else 0,
        metrics_table=render_metrics_table(metrics),
        filters_table=render_filters_table(metrics.get("filters", {})),
        cluster_table=cluster_size_table(adata),
        panels="\n".join(p for p in panels if p),
        generated_by="nf-scrnaseq report",
    )
    args.out.write_text(html)
    print(f"[make_report] wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

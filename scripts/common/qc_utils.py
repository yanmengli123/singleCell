# -*- coding: utf-8 -*-
"""QC utilities for single-cell analysis."""
import numpy as np


def calculate_qc_metrics(adata):
    """Calculate per-cell QC metrics."""
    # Number of genes per cell
    adata.obs['n_genes_by_counts'] = np.asarray(adata.X.sum(axis=1)).flatten()

    # Total counts per cell
    adata.obs['total_counts'] = np.asarray(adata.X.sum(axis=1)).flatten()

    return adata


def calculate_pct_mt(adata):
    """Calculate mitochondrial percentage for each cell."""
    if 'mt' not in adata.var.columns:
        adata.var['mt'] = adata.var.index.str.startswith('MT-')

    mt_genes = adata.var['mt']
    if mt_genes.sum() == 0:
        adata.obs['pct_counts_mt'] = 0.0
    else:
        adata.obs['pct_counts_mt'] = (
            adata[:, mt_genes].X.sum(axis=1).A1 / adata.obs['total_counts'] * 100
        )
    return adata


def filter_cells(adata, min_genes=200, max_genes=2500, max_pct_mt=5):
    """Filter cells based on QC criteria."""
    n_before = adata.n_obs

    adata = adata[
        (adata.obs['n_genes_by_counts'] >= min_genes) &
        (adata.obs['n_genes_by_counts'] <= max_genes) &
        (adata.obs['pct_counts_mt'] <= max_pct_mt)
    ].copy()

    n_after = adata.n_obs
    print(f"  Cells filtered: {n_before} -> {n_after} (removed {n_before - n_after})")
    return adata


def filter_genes(adata, min_cells=3):
    """Filter genes based on cell count."""
    n_before = adata.n_vars
    import scanpy as sc
    sc.pp.filter_genes(adata, min_cells=min_cells)
    n_after = adata.n_vars
    print(f"  Genes filtered: {n_before} -> {n_after} (removed {n_before - n_after})")
    return adata
# -*- coding: utf-8 -*-
"""QC utilities for single-cell analysis."""
import numpy as np
import scipy.sparse as sp


def calculate_qc_metrics(adata):
    """Calculate per-cell QC metrics.

    - n_genes_by_counts: number of genes with >0 counts per cell
    - total_counts: sum of all counts per cell
    """
    X = adata.X

    if sp.issparse(X):
        # For sparse matrix: sum axis=1 gives total counts
        # (X > 0).sum(axis=1) counts genes with >0
        n_genes_by_counts = np.array((X > 0).sum(axis=1)).flatten()
        total_counts = np.array(X.sum(axis=1)).flatten()
    else:
        # For dense array/matrix
        n_genes_by_counts = np.sum(X > 0, axis=1)
        total_counts = np.sum(X, axis=1)

    adata.obs['n_genes_by_counts'] = n_genes_by_counts
    adata.obs['total_counts'] = total_counts

    return adata


def calculate_pct_mt(adata):
    """Calculate mitochondrial percentage for each cell."""
    if 'mt' not in adata.var.columns:
        adata.var['mt'] = adata.var.index.str.startswith('MT-')

    mt_genes = adata.var['mt']
    if mt_genes.sum() == 0:
        adata.obs['pct_counts_mt'] = 0.0
    else:
        X_mt = adata[:, mt_genes].X
        if sp.issparse(X_mt):
            mt_counts = np.array(X_mt.sum(axis=1)).flatten()
        else:
            mt_counts = np.sum(X_mt, axis=1)
        adata.obs['pct_counts_mt'] = mt_counts / adata.obs['total_counts'] * 100

    return adata
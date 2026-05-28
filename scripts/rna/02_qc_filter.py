#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 2: QC Filtering
Adds QC metrics to obs/var, filters cells and genes, saves processed data.
"""
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path

# QC thresholds
MIN_GENES_PER_CELL = 200
MAX_GENES_PER_CELL = 2500
MAX_PCT_MT = 5
MIN_CELLS_PER_GENE = 3

def main():
    input_file = Path("data/raw/pbmc3k_raw.h5ad")
    output_file = Path("data/processed/pbmc3k_qc.h5ad")
    qc_summary_file = Path("results/tables/qc_summary.csv")

    print("=" * 50)
    print("STEP 2: QC FILTERING")
    print("=" * 50)

    # Load data
    adata = sc.read_h5ad(input_file)
    print(f"\nInput: {adata.n_obs} cells x {adata.n_vars} genes")

    # === Calculate QC metrics ===
    print("\n[1] Calculating QC metrics...")

    # Cell-level metrics
    adata.obs['n_genes_by_counts'] = np.asarray(adata.X.sum(axis=1)).flatten()
    adata.obs['total_counts'] = np.asarray(adata.X.sum(axis=1)).flatten()

    # Mitochondrial genes (MT- prefix)
    adata.var['mt'] = adata.var.index.str.startswith('MT-')
    adata.obs['pct_counts_mt'] = (
        adata[:, adata.var['mt']].X.sum(axis=1).A1 / adata.obs['total_counts'] * 100
    )

    print(f"   Cells with high mitochondrial %: {(adata.obs['pct_counts_mt'] > MAX_PCT_MT).sum()}")
    print(f"   Cells with too few genes: {(adata.obs['n_genes_by_counts'] < MIN_GENES_PER_CELL).sum()}")
    print(f"   Cells with too many genes: {(adata.obs['n_genes_by_counts'] > MAX_GENES_PER_CELL).sum()}")

    # === Filter cells ===
    print("\n[2] Filtering cells...")
    before_cells = adata.n_obs
    adata = adata[
        (adata.obs['n_genes_by_counts'] >= MIN_GENES_PER_CELL) &
        (adata.obs['n_genes_by_counts'] <= MAX_GENES_PER_CELL) &
        (adata.obs['pct_counts_mt'] <= MAX_PCT_MT)
    ].copy()
    after_cells = adata.n_obs
    print(f"   Cells before: {before_cells}, after: {after_cells} (removed {before_cells - after_cells})")

    # === Filter genes ===
    print("\n[3] Filtering genes...")
    before_genes = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
    after_genes = adata.n_vars
    print(f"   Genes before: {before_genes}, after: {after_genes} (removed {before_genes - after_genes})")

    # === Save raw counts in layers ===
    print("\n[4] Saving raw counts to layers...")
    adata.layers['counts'] = adata.X.copy()

    # === Save QC summary ===
    print("\n[5] Saving QC summary...")
    qc_summary_file.parent.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame({
        'metric': ['original_cells', 'filtered_cells', 'original_genes',
                   'filtered_genes', 'min_genes_per_cell', 'max_genes_per_cell',
                   'max_pct_mt', 'min_cells_per_gene'],
        'value': [before_cells, after_cells, before_genes, after_genes,
                  MIN_GENES_PER_CELL, MAX_GENES_PER_CELL, MAX_PCT_MT, MIN_CELLS_PER_GENE]
    })
    summary.to_csv(qc_summary_file, index=False)
    print(f"   Saved to: {qc_summary_file}")

    # === Save processed data ===
    print("\n[6] Saving filtered data...")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_file)
    print(f"   Saved to: {output_file}")

    print(f"\nFinal: {adata.n_obs} cells x {adata.n_vars} genes")
    print("\n" + "=" * 50)
    print("STEP 2 COMPLETE")
    print("=" * 50)

if __name__ == "__main__":
    main()
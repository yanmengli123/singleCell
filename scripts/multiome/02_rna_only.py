#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 2: Multiome RNA-only analysis
"""
import scanpy as sc
import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, str(Path(__file__).parent.parent))
from common.qc_utils import calculate_qc_metrics, calculate_pct_mt

PUBLICATION_DPI = 300

def main():
    input_file = Path("data/raw/multiome/rna.h5ad")
    output_file = Path("data/processed/multiome/rna_processed.h5ad")
    output_dir = Path("results/figures/multiome")

    print("=" * 50)
    print("STEP 2: MULTIOME RNA-ONLY ANALYSIS")
    print("=" * 50)

    adata = sc.read_h5ad(input_file)
    print(f"\nInput: {adata.n_obs} cells x {adata.n_vars} genes")

    # QC
    print("\n[1] QC...")
    calculate_qc_metrics(adata)
    calculate_pct_mt(adata)

    print(f"  n_genes_by_counts median: {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"  total_counts median: {adata.obs['total_counts'].median():.0f}")

    # Filter
    print("\n[2] Filtering cells and genes...")
    before = adata.n_obs
    adata = adata[
        (adata.obs['n_genes_by_counts'] >= 200) &
        (adata.obs['n_genes_by_counts'] <= 5000) &
        (adata.obs['pct_counts_mt'] <= 20)
    ].copy()
    print(f"  Cells: {before} -> {adata.n_obs}")

    sc.pp.filter_genes(adata, min_cells=3)
    print(f"  Genes: {adata.n_vars}")

    # Normalize
    print("\n[3] Normalizing...")
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)

    # HVG
    print("\n[4] Highly Variable Genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    print(f"  HVG: {adata.var['highly_variable'].sum()}")

    # Scale, PCA
    print("\n[5] PCA...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30, random_state=0)

    # Neighbors, UMAP, Clustering
    print("\n[6] Clustering...")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30, random_state=0)
    sc.tl.umap(adata, random_state=0)
    sc.tl.leiden(adata, resolution=0.5, random_state=0, flavor='igraph')
    print(f"  Clusters: {adata.obs['leiden'].nunique()}")

    # Save
    print("\n[7] Saving...")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_file)
    print(f"  Saved to: {output_file}")

    # Plot
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color='leiden', ax=ax, show=False, frameon=False,
               title='Multiome RNA-only UMAP', legend_loc='on data')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    fig.savefig(output_dir / 'rna_only_umap.png', dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_dir / 'rna_only_umap.png'}")

    print("\n" + "=" * 50)
    print("STEP 2 COMPLETE")
    print("=" * 50)

if __name__ == "__main__":
    main()
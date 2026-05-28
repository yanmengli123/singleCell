#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 3: Preprocess, PCA, Neighbors, Leiden, UMAP
"""
import scanpy as sc
import anndata as ad
from pathlib import Path
import matplotlib.pyplot as plt

# Parameters
N_TOP_GENES = 2000
N_PCS = 30
NEIGHBORS_N = 10
LEIDEN_RESOLUTION = 0.5
RANDOM_STATE = 0

def main():
    input_file = Path("data/processed/pbmc3k_qc.h5ad")
    output_file = Path("data/processed/pbmc3k_processed.h5ad")

    print("=" * 50)
    print("STEP 3: PREPROCESS & CLUSTER")
    print("=" * 50)

    adata = sc.read_h5ad(input_file)
    print(f"\nInput: {adata.n_obs} cells x {adata.n_vars} genes")

    # === Normalize and log transform ===
    print("\n[1] Normalizing and log-transforming...")
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)

    # === Highly variable genes ===
    print("\n[2] Selecting highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES)
    n_hvg = adata.var['highly_variable'].sum()
    print(f"   HVG selected: {n_hvg}")

    # === Scale ===
    print("\n[3] Scaling data...")
    sc.pp.scale(adata, max_value=10)

    # === PCA ===
    print("\n[4] Running PCA...")
    sc.tl.pca(adata, n_comps=N_PCS, random_state=RANDOM_STATE)
    print(f"   PCA shape: {adata.obsm['X_pca'].shape}")

    # === Variance plot ===
    fig, ax = plt.subplots(1, 1, figsize=(4, 3))
    ax.plot(range(1, N_PCS + 1), adata.varm['PCs'].var(axis=0)[:N_PCS], 'o-')
    ax.set_xlabel('PC'); ax.set_ylabel('Variance')
    ax.set_title('PCA Variance')
    fig.tight_layout()
    fig.savefig("results/figures/pca_variance.png", dpi=100)
    plt.close()
    print("   Saved: results/figures/pca_variance.png")

    # === Neighbors ===
    print("\n[5] Computing neighbors...")
    sc.pp.neighbors(adata, n_neighbors=NEIGHBORS_N, n_pcs=N_PCS, random_state=RANDOM_STATE)

    # === UMAP ===
    print("\n[6] Running UMAP...")
    sc.tl.umap(adata, random_state=RANDOM_STATE)

    # === Leiden clustering ===
    print("\n[7] Leiden clustering...")
    sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, random_state=RANDOM_STATE)
    n_clusters = adata.obs['leiden'].nunique()
    print(f"   Clusters: {n_clusters}")

    # === Save ===
    print("\n[8] Saving results...")
    adata.write_h5ad(output_file)
    print(f"   Saved to: {output_file}")

    # === Cluster size distribution ===
    print("\n[9] Cluster distribution:")
    print(adata.obs['leiden'].value_counts().sort_index())

    # === UMAP plot ===
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    sc.pl.umap(adata, color='leiden', ax=ax, show=False)
    fig.savefig("results/figures/umap_leiden.png", dpi=100, bbox_inches='tight')
    plt.close()
    print("\n   Saved: results/figures/umap_leiden.png")

    print("\n" + "=" * 50)
    print("STEP 3 COMPLETE")
    print("=" * 50)

if __name__ == "__main__":
    main()
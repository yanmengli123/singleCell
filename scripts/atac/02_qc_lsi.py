#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 2: ATAC QC and LSI dimensionality reduction
"""
import scanpy as sc
import muon as mu
from muon import atac as atac
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

PUBLICATION_DPI = 300

def setup_publication_style():
    """Configure matplotlib for publication-quality figures."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10,
        'axes.linewidth': 1.5,
        'axes.titlesize': 11,
        'axes.labelsize': 10,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'xtick.major.width': 1.2,
        'ytick.major.width': 1.2,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.dpi': PUBLICATION_DPI,
        'savefig.bbox': 'tight',
        'savefig.facecolor': 'white',
        'legend.frameon': False,
    })


def main():
    input_file = Path("data/raw/atac/pbmc10k_atac_raw.h5ad")
    output_file = Path("data/processed/atac/pbmc10k_atac_qc.h5ad")
    output_dir = Path("results/figures/atac")

    print("=" * 50)
    print("STEP 2: ATAC QC & LSI")
    print("=" * 50)

    setup_publication_style()

    # Load raw ATAC data
    print("\n[1] Loading raw ATAC data...")
    adata = sc.read_h5ad(input_file)
    print(f"Input: {adata.n_obs} cells x {adata.n_vars} peaks")

    # Calculate QC metrics
    print("\n[2] Calculating QC metrics...")
    n_peaks = np.array(adata.X.sum(axis=1)).flatten()
    adata.obs['n_peaks'] = n_peaks
    adata.obs['total_counts'] = n_peaks

    # TSS enrichment (simplified: use peak counts distribution)
    # In real analysis, we'd use fragment file and TSS annotation
    print(f"  Median peaks per cell: {adata.obs['n_peaks'].median():.0f}")
    print(f"  Total cells: {adata.n_obs}")

    # QC filtering - similar to RNA but for ATAC
    print("\n[3] QC Filtering...")
    min_peaks = 1000
    max_peaks = 100000

    before = adata.n_obs
    adata = adata[
        (adata.obs['n_peaks'] >= min_peaks) &
        (adata.obs['n_peaks'] <= max_peaks)
    ].copy()
    after = adata.n_obs
    print(f"  Cells: {before} -> {after} (removed {before - after})")

    # Also filter low-quality cells using metadata if available
    if 'is__cell_barcode' in adata.obs.columns:
        n_before = adata.n_obs
        adata = adata[adata.obs['is__cell_barcode'] == 1].copy()
        print(f"  Cell barcode filter: {n_before} -> {adata.n_obs}")

    # Normalize and log transform
    print("\n[4] Normalizing and log-transforming...")
    adata.layers['counts'] = adata.X.copy()

    # TF-IDF normalization for ATAC
    # TF: term frequency (normalize per cell)
    # IDF: inverse document frequency (weight peaks by frequency)
    from sklearn.feature_extraction.text import TfidfTransformer
    tfidf = TfidfTransformer()
    X_tfidf = tfidf.fit_transform(adata.X.T).T  # transpose for TF-IDF

    adata.X = X_tfidf.tocsr()
    print(f"  TF-IDF normalized: {X_tfidf.shape}")

    # LSI (Latent Semantic Indexing) - SVD on TF-IDF
    print("\n[5] LSI (SVD) dimensionality reduction...")
    from sklearn.decomposition import TruncatedSVD
    n_comps = 30
    svd = TruncatedSVD(n_components=n_comps, random_state=0)
    X_lsi = svd.fit_transform(adata.X)

    adata.obsm['X_lsi'] = X_lsi
    adata.varm['LSI_components'] = svd.components_.T
    print(f"  LSI shape: {X_lsi.shape}")
    print(f"  Variance explained: {svd.explained_variance_ratio_.sum():.4%}")

    # Neighbors and UMAP
    print("\n[6] Computing neighbors and UMAP...")
    sc.pp.neighbors(adata, use_rep='X_lsi', n_neighbors=10, random_state=0)
    sc.tl.umap(adata, random_state=0)

    # Leiden clustering
    print("\n[7] Leiden clustering...")
    sc.tl.leiden(adata, resolution=0.5, random_state=0, flavor='igraph')
    print(f"  Clusters: {adata.obs['leiden'].nunique()}")

    # Save
    print("\n[8] Saving processed ATAC data...")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_file)
    print(f"  Saved to: {output_file}")

    # Publication figures
    print("\n[9] Generating publication figures...")
    output_dir.mkdir(parents=True, exist_ok=True)

    # QC violin
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    for ax, col, label, color in zip(axes,
                                      ['n_peaks', 'total_counts'],
                                      ['Peaks per cell', 'Total counts'],
                                      ['steelblue', 'lightcoral']):
        parts = ax.violinplot(adata.obs[col].values, showmeans=False, showmedians=True)
        for pc in parts['bodies']:
            pc.set_facecolor(color)
            pc.set_alpha(0.7)
        parts['cmedians'].set_color('black')
        ax.set_ylabel(label)
        ax.set_xticks([1])
        ax.set_xticklabels(['All cells'])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    plt.tight_layout()
    fig.savefig(output_dir / 'atac_qc_violin.png', dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()

    # Variance explained
    fig, ax = plt.subplots(figsize=(5, 4))
    var_ratio = svd.explained_variance_ratio_[:n_comps]
    cumulative = np.cumsum(var_ratio)
    x = np.arange(1, n_comps + 1)
    ax.bar(x, var_ratio, color='steelblue', alpha=0.7, label='Individual')
    ax.plot(x, cumulative, 'o-', color='darkblue', markersize=3, label='Cumulative')
    ax.set_xlabel('SVD Component')
    ax.set_ylabel('Variance Explained')
    ax.set_title('LSI Variance Explained')
    ax.legend(loc='center right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    fig.savefig(output_dir / 'lsi_variance.png', dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()

    # UMAP by Leiden
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.umap(adata, color='leiden', ax=ax, show=False, frameon=False,
               title='ATAC Leiden Clusters', legend_loc='on data')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    fig.savefig(output_dir / 'umap_leiden.png', dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()

    # Cluster distribution
    print("\n[10] Cluster distribution:")
    print(adata.obs['leiden'].value_counts().sort_index())

    print("\n" + "=" * 50)
    print("STEP 2 COMPLETE")
    print("=" * 50)


if __name__ == "__main__":
    main()
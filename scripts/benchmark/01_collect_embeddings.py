#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Benchmark: Compare different multi-omics integration methods
Compares RNA-only, ATAC-only, and simulated integration results
"""
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

PUBLICATION_DPI = 300


def setup_publication_style():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10,
        'axes.linewidth': 1.5,
        'axes.titlesize': 11,
        'axes.labelsize': 10,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.dpi': PUBLICATION_DPI,
        'savefig.bbox': 'tight',
    })


def compute_integration_metrics(adata, embedding_key='X_pca', cell_type_key='cell_type'):
    """Compute basic integration metrics."""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.model_selection import cross_val_score

    X = adata.obsm[embedding_key]
    y = adata.obs[cell_type_key].values

    # ARI using KNN
    knn = KNeighborsClassifier(n_neighbors=10)
    scores = cross_val_score(knn, X, y, cv=5, scoring='accuracy')

    # Silhouette score
    from sklearn.metrics import silhouette_score
    labels = adata.obs['leiden'].values
    sil = silhouette_score(X[:1000], labels[:1000]) if len(X) > 1000 else silhouette_score(X, labels)

    return {
        'knn_accuracy': scores.mean(),
        'silhouette_score': sil
    }


def main():
    setup_publication_style()

    print("=" * 50)
    print("BENCHMARK: COMPARISON OF METHODS")
    print("=" * 50)

    output_dir = Path("results/figures/benchmark")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load GLUE tutorial data
    print("\n[1] Loading GLUE tutorial data...")
    adata_rna = sc.read_h5ad("data/processed/glue/rna_prepared.h5ad")
    adata_atac = sc.read_h5ad("data/processed/glue/atac_prepared.h5ad")

    print(f"RNA: {adata_rna.n_obs} cells x {adata_rna.n_vars} genes")
    print(f"ATAC: {adata_atac.n_obs} cells x {adata_atac.n_vars} peaks")

    # Method 1: RNA-only baseline
    print("\n[2] RNA-only baseline...")
    sc.pp.normalize_total(adata_rna, target_sum=10000)
    sc.pp.log1p(adata_rna)
    sc.pp.highly_variable_genes(adata_rna, n_top_genes=2000)
    adata_rna = adata_rna[:, adata_rna.var['highly_variable']].copy()
    sc.pp.scale(adata_rna, max_value=10)
    sc.tl.pca(adata_rna, n_comps=30, random_state=0)
    sc.pp.neighbors(adata_rna, n_neighbors=10, random_state=0)
    sc.tl.umap(adata_rna, random_state=0)
    sc.tl.leiden(adata_rna, resolution=0.5, random_state=0, flavor='igraph')
    print(f"  RNA clusters: {adata_rna.obs['leiden'].nunique()}")

    # Method 2: ATAC-only baseline (using LSI)
    print("\n[3] ATAC-only baseline (LSI)...")
    from sklearn.feature_extraction.text import TfidfTransformer
    from sklearn.decomposition import TruncatedSVD

    # TF-IDF + LSI
    tfidf = TfidfTransformer()
    X_tfidf = tfidf.fit_transform(adata_atac.X.T).T
    svd = TruncatedSVD(n_components=30, random_state=0)
    X_lsi = svd.fit_transform(X_tfidf.tocsr())
    adata_atac.obsm['X_lsi'] = X_lsi

    sc.pp.neighbors(adata_atac, use_rep='X_lsi', n_neighbors=10, random_state=0)
    sc.tl.umap(adata_atac, random_state=0)
    sc.tl.leiden(adata_atac, resolution=0.5, random_state=0, flavor='igraph')
    print(f"  ATAC clusters: {adata_atac.obs['leiden'].nunique()}")

    # Compute metrics
    print("\n[4] Computing metrics...")
    metrics = {}

    if 'cell_type' in adata_rna.obs.columns:
        metrics['RNA'] = compute_integration_metrics(adata_rna, 'X_pca', 'cell_type')
        print(f"  RNA KNN accuracy: {metrics['RNA']['knn_accuracy']:.4f}")

    if 'cell_type' in adata_atac.obs.columns:
        metrics['ATAC'] = compute_integration_metrics(adata_atac, 'X_lsi', 'cell_type')
        print(f"  ATAC KNN accuracy: {metrics['ATAC']['knn_accuracy']:.4f}")

    # Generate publication figures
    print("\n[5] Generating publication figures...")

    # Figure 1: RNA-only UMAP
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    sc.pl.umap(adata_rna, color='cell_type', ax=ax, show=False, frameon=False,
               title='RNA-only UMAP', legend_loc='right margin', legend_fontsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax = axes[1]
    sc.pl.umap(adata_atac, color='cell_type', ax=ax, show=False, frameon=False,
               title='ATAC-only UMAP', legend_loc='right margin', legend_fontsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_dir / 'rna_vs_atac_umap.png', dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_dir / 'rna_vs_atac_umap.png'}")

    # Figure 2: Metrics comparison
    if metrics:
        fig, ax = plt.subplots(figsize=(6, 4))
        methods = list(metrics.keys())
        acc_values = [metrics[m].get('knn_accuracy', 0) for m in methods]

        bars = ax.bar(methods, acc_values, color=['steelblue', 'lightcoral'])
        ax.set_ylabel('KNN Accuracy')
        ax.set_title('Integration Benchmark: Cell Type Prediction')
        ax.set_ylim(0, 1.0)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Add value labels on bars
        for bar, val in zip(bars, acc_values):
            ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.02,
                    f'{val:.3f}', ha='center', va='bottom', fontsize=9)

        plt.tight_layout()
        fig.savefig(output_dir / 'benchmark_metrics.png', dpi=PUBLICATION_DPI, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_dir / 'benchmark_metrics.png'}")

    # Save metrics
    print("\n[6] Saving benchmark results...")
    metrics_df = pd.DataFrame(metrics).T
    metrics_df.to_csv("results/tables/benchmark/benchmark_summary.csv")
    print("  Saved: results/tables/benchmark/benchmark_summary.csv")

    print("\n" + "=" * 50)
    print("BENCHMARK COMPLETE")
    print("=" * 50)


if __name__ == "__main__":
    main()
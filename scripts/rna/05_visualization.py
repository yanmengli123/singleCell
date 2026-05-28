#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 5: Publication-quality visualization for RNA pipeline
"""
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Publication style settings
PUBLICATION_DPI = 300
FONT_SIZE = 10
AXIS_LINEWIDTH = 1.5
TICK_WIDTH = 1.2

def setup_publication_style():
    """Configure matplotlib for publication-quality figures."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': FONT_SIZE,
        'axes.linewidth': AXIS_LINEWIDTH,
        'axes.titlesize': FONT_SIZE + 1,
        'axes.labelsize': FONT_SIZE,
        'xtick.labelsize': FONT_SIZE - 1,
        'ytick.labelsize': FONT_SIZE - 1,
        'xtick.major.width': TICK_WIDTH,
        'ytick.major.width': TICK_WIDTH,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.dpi': PUBLICATION_DPI,
        'savefig.bbox': 'tight',
        'savefig.facecolor': 'white',
        'legend.frameon': False,
        'legend.fontsize': FONT_SIZE - 1,
    })


def plot_qc_violin(adata, output_path):
    """Publication QC violin plots."""
    setup_publication_style()
    fig, axes = plt.subplots(1, 3, figsize=(10, 4))

    metrics = [
        ('n_genes_by_counts', 'Genes per cell', 'lightblue'),
        ('total_counts', 'Total UMI counts', 'lightcoral'),
        ('pct_counts_mt', 'Mitochondrial %', 'lightgreen')
    ]

    for ax, (col, label, color) in zip(axes, metrics):
        if col in adata.obs.columns:
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
    fig.savefig(output_path, dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_qc_scatter(adata, output_path):
    """Publication QC scatter plots."""
    setup_publication_style()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Total counts vs n_genes
    ax = axes[0]
    ax.scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'],
               c='steelblue', alpha=0.5, s=5, rasterized=True)
    ax.set_xlabel('Total UMI counts')
    ax.set_ylabel('Genes per cell')
    ax.set_title('QC: Counts vs Genes')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # pct_mt vs total_counts
    ax = axes[1]
    ax.scatter(adata.obs['total_counts'], adata.obs['pct_counts_mt'],
               c='tomato', alpha=0.5, s=5, rasterized=True)
    ax.set_xlabel('Total UMI counts')
    ax.set_ylabel('Mitochondrial %')
    ax.set_title('QC: Counts vs MT%')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_pca_variance(adata, n_pcs, output_path):
    """Publication PCA variance explained plot."""
    setup_publication_style()
    fig, ax = plt.subplots(figsize=(5, 4))

    variance_ratio = adata.uns['pca']['variance_ratio'][:n_pcs]
    cumulative = np.cumsum(variance_ratio)

    x = np.arange(1, n_pcs + 1)
    ax.bar(x, variance_ratio, color='steelblue', alpha=0.7, label='Individual')
    ax.plot(x, cumulative, 'o-', color='darkblue', markersize=3, label='Cumulative')
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Variance Explained')
    ax.set_title('PCA Variance Explained')
    ax.legend(loc='center right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0.5, n_pcs + 0.5)

    plt.tight_layout()
    fig.savefig(output_path, dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_umap_leiden(adata, output_path):
    """Publication UMAP with Leiden clusters."""
    setup_publication_style()
    fig, ax = plt.subplots(figsize=(7, 6))

    n_clusters = adata.obs['leiden'].nunique()
    palette = sc.pl.palettes.default_102 if n_clusters > 20 else None

    sc.pl.umap(adata, color='leiden', ax=ax, show=False, frameon=False,
               title='Leiden Clusters', legend_loc='on data' if n_clusters <= 10 else 'right margin',
               palette=palette)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_umap_celltype(adata, output_path):
    """Publication UMAP with cell types."""
    setup_publication_style()
    fig, ax = plt.subplots(figsize=(8, 6))

    cell_type_col = 'cell_type' if 'cell_type' in adata.obs.columns else 'leiden'
    sc.pl.umap(adata, color=cell_type_col, ax=ax, show=False, frameon=False,
               title='Cell Type Annotation', legend_loc='right margin')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_path, dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_marker_dotplot(adata, marker_dict, groupby, output_path):
    """Publication marker gene dotplot."""
    setup_publication_style()
    fig = plt.figure(figsize=(12, 6))

    available_markers = {}
    for ct, genes in marker_dict.items():
        avail = [g for g in genes if g in adata.var_names]
        if avail:
            available_markers[ct] = avail

    if available_markers:
        sc.pl.dotplot(adata, available_markers, groupby=groupby, show=False,
                      standard_scale='var', return_fig=True)
        plt.tight_layout()
        fig.savefig(output_path, dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_marker_heatmap(adata, marker_dict, groupby, output_path):
    """Publication marker gene heatmap."""
    setup_publication_style()

    available_markers = []
    for ct, genes in marker_dict.items():
        avail = [g for g in genes if g in adata.var_names]
        available_markers.extend(avail)

    if available_markers:
        fig = sc.pl.heatmap(adata, available_markers, groupby=groupby,
                            show=False, swap_axes=True, figsize=(10, 8),
                            show_gene_labels=True)
        if hasattr(fig, 'savefig'):
            fig.savefig(output_path, dpi=PUBLICATION_DPI, bbox_inches='tight')
        else:
            plt.gcf().savefig(output_path, dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_rank_genes(adata, output_path):
    """Publication ranked genes per cluster."""
    setup_publication_style()
    result = sc.pl.rank_genes_groups(adata, n_genes=10, show=False, return_fig=True,
                                     sharey=False, fontsize=9)
    plt.tight_layout()
    # result could be Axes or list of Axes
    if isinstance(result, (list, np.ndarray)):
        fig = result[0].get_figure() if hasattr(result[0], 'get_figure') else plt.gcf()
    elif hasattr(result, 'get_figure'):
        fig = result.get_figure()
    else:
        fig = plt.gcf()
    fig.savefig(output_path, dpi=PUBLICATION_DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    print("=" * 50)
    print("STEP 5: PUBLICATION VISUALIZATION")
    print("=" * 50)

    input_file = Path("data/processed/rna/pbmc3k_annotated.h5ad")
    output_dir = Path("results/figures/rna")
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(input_file)
    print(f"\nInput: {adata.n_obs} cells x {adata.n_vars} genes")

    # PBMC marker genes
    MARKER_GENES = {
        'CD4 T': ['CD3D', 'CD3E', 'IL7R', 'CCR7'],
        'CD8 T': ['CD8A', 'CD8B'],
        'B cells': ['MS4A1', 'CD79A'],
        'NK cells': ['NKG7', 'GNLY'],
        'Monocytes': ['S100A8', 'S100A9', 'FCGR3A', 'LST1'],
        'DC': ['FCER1A', 'CST3'],
        'Platelets': ['PPBP']
    }

    # Generate all publication figures
    print("\n[1] QC Violin Plot...")
    plot_qc_violin(adata, output_dir / 'qc_violin.png')

    print("\n[2] QC Scatter Plot...")
    plot_qc_scatter(adata, output_dir / 'qc_scatter.png')

    print("\n[3] PCA Variance...")
    plot_pca_variance(adata, 30, output_dir / 'pca_variance.png')

    print("\n[4] UMAP Leiden...")
    plot_umap_leiden(adata, output_dir / 'umap_leiden.png')

    print("\n[5] UMAP Cell Type...")
    plot_umap_celltype(adata, output_dir / 'umap_celltype.png')

    print("\n[6] Marker Dotplot...")
    plot_marker_dotplot(adata, MARKER_GENES, 'leiden', output_dir / 'marker_dotplot.png')

    print("\n[7] Marker Heatmap...")
    plot_marker_heatmap(adata, MARKER_GENES, 'leiden', output_dir / 'marker_heatmap.png')

    print("\n[8] Rank Genes Groups...")
    if 'rank_genes_groups' in adata.uns:
        plot_rank_genes(adata, output_dir / 'rank_genes.png')

    print("\n" + "=" * 50)
    print("STEP 5 COMPLETE")
    print("=" * 50)

if __name__ == "__main__":
    main()
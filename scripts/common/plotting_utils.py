# -*- coding: utf-8 -*-
"""Plotting utilities for single-cell analysis."""
import matplotlib.pyplot as plt
import scanpy as sc
from pathlib import Path


def setup_figure_style():
    """Set up matplotlib style for single-cell plots."""
    sc.settings.set_figure_params(dpi=80, facecolor='white', frameon=False)


def save_figure(fig, path, dpi=100):
    """Save figure to path, creating directory if needed."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def plot_qc_violin(adata, qc_metrics, path):
    """Plot QC metric violin plots."""
    fig = sc.pl.violin(adata, qc_metrics, show=False, return_fig=True)
    save_figure(fig, path)


def plot_qc_scatter(adata, x_metric, y_metric, color_metric, path):
    """Plot scatter plot for QC metrics."""
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    ax.scatter(adata.obs[x_metric], adata.obs[y_metric],
               c=adata.obs[color_metric], alpha=0.5, s=5)
    ax.set_xlabel(x_metric)
    ax.set_ylabel(y_metric)
    fig.tight_layout()
    save_figure(fig, path)


def plot_umap(adata, color, title, path):
    """Plot UMAP colored by given column."""
    fig = sc.pl.umap(adata, color=color, title=title, show=False, return_fig=True)
    save_figure(fig, path)
# -*- coding: utf-8 -*-
"""Tests for data integrity checks."""
import os
import sys
from pathlib import Path
import numpy as np


def test_raw_h5ad_exists():
    """Test that raw h5ad file exists at correct path."""
    path = Path("data/raw/rna/pbmc3k_raw.h5ad")
    assert path.exists(), f"Raw h5ad not found: {path}"


def test_paths_consistent_with_config():
    """Test that script paths match config paths."""
    import yaml
    with open("configs/01_pbmc3k_rna.yaml", 'r') as f:
        config = yaml.safe_load(f)
    input_path = config['paths']['input']
    assert Path(input_path).exists(), f"Config input path not found: {input_path}"


def test_processed_h5ad_files_exist():
    """Test that all processed h5ad files exist after running pipeline."""
    files = [
        "data/processed/rna/pbmc3k_qc.h5ad",
        "data/processed/rna/pbmc3k_processed.h5ad",
        "data/processed/rna/pbmc3k_annotated.h5ad",
    ]
    for f in files:
        assert Path(f).exists(), f"Processed h5ad not found: {f}"


def test_qc_calculation_correct():
    """Test that n_genes_by_counts != total_counts after QC fix."""
    import scanpy as sc
    adata = sc.read_h5ad("data/processed/rna/pbmc3k_qc.h5ad")

    assert 'n_genes_by_counts' in adata.obs, "n_genes_by_counts not in obs"
    assert 'total_counts' in adata.obs, "total_counts not in obs"

    # They should NOT be equal for most cells
    n_equal = (adata.obs['n_genes_by_counts'] == adata.obs['total_counts']).sum()
    pct_equal = n_equal / adata.n_obs * 100
    assert pct_equal < 50, f"n_genes_by_counts == total_counts for {pct_equal:.1f}% of cells - QC may be wrong"

    # n_genes_by_counts should be much smaller than total_counts
    median_ngenes = adata.obs['n_genes_by_counts'].median()
    median_counts = adata.obs['total_counts'].median()
    assert median_ngenes < median_counts / 2, \
        f"n_genes_by_counts ({median_ngenes:.0f}) too close to total_counts ({median_counts:.0f})"


def test_pct_counts_mt_valid():
    """Test that pct_counts_mt is in valid range 0-100."""
    import scanpy as sc
    adata = sc.read_h5ad("data/processed/rna/pbmc3k_qc.h5ad")

    if 'pct_counts_mt' in adata.obs:
        assert adata.obs['pct_counts_mt'].min() >= 0, "pct_counts_mt below 0"
        assert adata.obs['pct_counts_mt'].max() <= 100, "pct_counts_mt above 100"


def test_counts_layer_exists():
    """Test that counts layer exists and is raw data."""
    import scanpy as sc
    adata = sc.read_h5ad("data/processed/rna/pbmc3k_qc.h5ad")

    assert 'counts' in adata.layers, "counts layer not found"
    # counts should be raw (non-negative integers)
    X_counts = adata.layers['counts']
    if hasattr(X_counts, 'toarray'):
        X_counts = X_counts.toarray()
    assert (X_counts >= 0).all(), "counts layer has negative values"
    assert (X_counts == X_counts.astype(int)).all(), "counts layer has non-integer values"


def test_pca_umap_leiden_exist():
    """Test that PCA, UMAP, and Leiden exist in processed data."""
    import scanpy as sc
    adata = sc.read_h5ad("data/processed/rna/pbmc3k_processed.h5ad")

    errors = []
    if 'X_pca' not in adata.obsm:
        errors.append("X_pca not in obsm")
    if 'X_umap' not in adata.obsm:
        errors.append("X_umap not in obsm")
    if 'leiden' not in adata.obs:
        errors.append("leiden not in obs")

    assert len(errors) == 0, f"Missing in processed: {', '.join(errors)}"


def test_cell_type_no_na():
    """Test that cell_type has no NA values."""
    import scanpy as sc
    adata = sc.read_h5ad("data/processed/rna/pbmc3k_annotated.h5ad")

    assert 'cell_type' in adata.obs, "cell_type not in obs"
    na_count = adata.obs['cell_type'].isna().sum()
    assert na_count == 0, f"cell_type has {na_count} NA values"

    tbd_count = (adata.obs['cell_type'] == 'TBD').sum()
    if tbd_count > 0:
        print(f"  WARNING: {tbd_count} cells have TBD annotation")


def test_cluster_annotation_covers_all_clusters():
    """Test that all leiden clusters have annotation in config."""
    import scanpy as sc
    import yaml
    adata = sc.read_h5ad("data/processed/rna/pbmc3k_processed.h5ad")

    with open("configs/pbmc3k_cluster_annotation.yaml", 'r') as f:
        config = yaml.safe_load(f)
    cluster_annot = config['cluster_annotation']

    leiden_clusters = set(adata.obs['leiden'].unique())
    annotated_clusters = set(cluster_annot.keys())

    missing = leiden_clusters - annotated_clusters
    assert len(missing) == 0, f"Clusters {missing} not in annotation config"


def test_readme_commands_exist():
    """Test that README commands reference existing files."""
    readme = Path("README.md").read_text()

    # Check scripts exist
    scripts = [
        "scripts/rna/00_check_environment.py",
        "scripts/rna/01_inspect_anndata.py",
        "scripts/rna/02_qc_filter.py",
        "scripts/rna/03_preprocess_cluster.py",
        "scripts/rna/04_marker_annotation.py",
    ]
    for s in scripts:
        assert Path(s).exists(), f"README references missing script: {s}"


def test_qc_summary_exists():
    """Test that QC summary CSV exists."""
    path = Path("results/tables/rna/qc_summary.csv")
    assert path.exists(), f"QC summary not found: {path}"


def test_marker_genes_exists():
    """Test that marker genes CSV exists."""
    path = Path("results/tables/rna/marker_genes.csv")
    assert path.exists(), f"Marker genes CSV not found: {path}"


if __name__ == "__main__":
    print("Running data integrity tests...\n")

    tests = [
        test_raw_h5ad_exists,
        test_paths_consistent_with_config,
        test_processed_h5ad_files_exist,
        test_pct_counts_mt_valid,
        test_counts_layer_exists,
        test_pca_umap_leiden_exist,
        test_cell_type_no_na,
        test_cluster_annotation_covers_all_clusters,
        test_readme_commands_exist,
        test_qc_summary_exists,
        test_marker_genes_exists,
    ]

    passed = 0
    failed = 0

    for test in tests:
        try:
            test()
            print(f"  PASS: {test.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {test.__name__} - {e}")
            failed += 1
        except Exception as e:
            print(f"  ERROR: {test.__name__} - {e}")
            failed += 1

    print(f"\n{passed} passed, {failed} failed")

    # Run QC-specific test separately (requires loading h5ad)
    print("\n--- Running QC calculation test (requires h5ad) ---")
    try:
        test_qc_calculation_correct()
        print("  PASS: test_qc_calculation_correct")
    except AssertionError as e:
        print(f"  FAIL: test_qc_calculation_correct - {e}")
        failed += 1
    except Exception as e:
        print(f"  ERROR: test_qc_calculation_correct - {e}")
        failed += 1
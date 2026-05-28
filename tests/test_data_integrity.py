# -*- coding: utf-8 -*-
"""Tests for data integrity checks."""
import os
import sys
from pathlib import Path


def test_raw_h5ad_exists():
    """Test that raw h5ad file exists."""
    path = Path("data/raw/rna/pbmc3k_raw.h5ad")
    assert path.exists(), f"Raw h5ad not found: {path}"


def test_processed_h5ad_files_exist():
    """Test that all processed h5ad files exist."""
    files = [
        "data/processed/rna/pbmc3k_qc.h5ad",
        "data/processed/rna/pbmc3k_processed.h5ad",
        "data/processed/rna/pbmc3k_annotated.h5ad",
    ]
    for f in files:
        assert Path(f).exists(), f"Processed h5ad not found: {f}"


def test_annotated_h5ad_has_required_fields():
    """Test that annotated h5ad has required fields."""
    import scanpy as sc
    adata = sc.read_h5ad("data/processed/rna/pbmc3k_annotated.h5ad")

    errors = []

    if adata.n_obs == 0:
        errors.append("0 cells in annotated h5ad")
    if adata.n_vars == 0:
        errors.append("0 genes in annotated h5ad")
    if 'X_pca' not in adata.obsm:
        errors.append("X_pca not in obsm")
    if 'X_umap' not in adata.obsm:
        errors.append("X_umap not in obsm")
    if 'leiden' not in adata.obs:
        errors.append("leiden not in obs")
    if 'cell_type' not in adata.obs:
        errors.append("cell_type not in obs")

    assert len(errors) == 0, f"Errors: {', '.join(errors)}"


def test_qc_summary_exists():
    """Test that QC summary CSV exists."""
    path = Path("results/tables/rna/qc_summary.csv")
    assert path.exists(), f"QC summary not found: {path}"


def test_marker_genes_exists():
    """Test that marker genes CSV exists."""
    path = Path("results/tables/rna/marker_genes.csv")
    assert path.exists(), f"Marker genes CSV not found: {path}"


if __name__ == "__main__":
    print("Running data integrity tests...")

    tests = [
        test_raw_h5ad_exists,
        test_processed_h5ad_files_exist,
        test_qc_summary_exists,
        test_marker_genes_exists,
    ]

    for test in tests:
        try:
            test()
            print(f"  PASS: {test.__name__}")
        except AssertionError as e:
            print(f"  FAIL: {test.__name__} - {e}")
        except Exception as e:
            print(f"  ERROR: {test.__name__} - {e}")

    print("\nNote: Run test_annotated_h5ad_has_required_fields separately with Python")
    print("as it requires loading the h5ad file.")
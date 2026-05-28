# -*- coding: utf-8 -*-
"""Validation utilities for single-cell analysis pipeline."""
import os


def validate_anndata(adata, name=""):
    """Validate AnnData object has required properties."""
    errors = []

    if adata.n_obs == 0:
        errors.append(f"{name}: 0 cells")
    if adata.n_vars == 0:
        errors.append(f"{name}: 0 genes")
    if adata.X is None:
        errors.append(f"{name}: X is None")

    return errors


def validate_file_exists(path, description=""):
    """Validate a file exists."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Required {description} not found: {path}")
    return True


def validate_h5ad(path, key=""):
    """Validate h5ad file can be loaded."""
    import scanpy as sc

    validate_file_exists(path, f"h5ad file {key}")

    try:
        adata = sc.read_h5ad(path)
        errors = validate_anndata(adata, key)
        if errors:
            raise ValueError(f"Invalid h5ad {key}: {', '.join(errors)}")
        return adata
    except Exception as e:
        raise RuntimeError(f"Failed to load h5ad {path}: {e}")


def validate_required_fields(adata, fields, location="obs"):
    """Validate required fields exist in adata."""
    missing = []
    if location == "obs":
        for f in fields:
            if f not in adata.obs.columns:
                missing.append(f"obs.{f}")
    elif location == "var":
        for f in fields:
            if f not in adata.var.columns:
                missing.append(f"var.{f}")
    elif location == "obsm":
        for f in fields:
            if f not in adata.obsm:
                missing.append(f"obsm[{f}]")
    elif location == "layers":
        for f in fields:
            if f not in adata.layers:
                missing.append(f"layers[{f}]")

    if missing:
        raise ValueError(f"Missing required fields: {', '.join(missing)}")

    return True
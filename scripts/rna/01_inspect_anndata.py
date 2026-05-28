#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 1: Inspect AnnData structure
Reads the raw h5ad and prints key information.
"""
import scanpy as sc
from pathlib import Path

def main():
    data_file = Path("data/raw/rna/pbmc3k_raw.h5ad")
    if not data_file.exists():
        print(f"ERROR: {data_file} not found!")
        return

    print("=" * 50)
    print("STEP 1: INSPECT AnnData")
    print("=" * 50)

    # Load data
    adata = sc.read_h5ad(data_file)

    # Basic shape
    print(f"\n1. Raw matrix shape: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"   Shape tuple: {adata.shape}")

    # Expression matrix
    print(f"\n2. Expression matrix (adata.X):")
    print(f"   Type: {type(adata.X)}")
    print(f"   Density: {adata.X.nnz / (adata.n_obs * adata.n_vars):.4%}")

    # obs (cell metadata)
    print(f"\n3. Cell metadata (adata.obs):")
    print(f"   Columns: {list(adata.obs.columns)}")
    print(f"   Shape: {adata.obs.shape}")
    if adata.obs.shape[1] > 0:
        print(adata.obs.head())

    # var (gene metadata)
    print(f"\n4. Gene metadata (adata.var):")
    print(f"   Columns: {list(adata.var.columns)}")
    print(f"   Shape: {adata.var.shape}")
    print(adata.var.head())

    # obsm (multi-dimensional matrices)
    print(f"\n5. Multi-dimensional matrices (adata.obsm):")
    print(f"   Keys: {list(adata.obsm.keys())}")

    # layers
    print(f"\n6. Layers (adata.layers):")
    print(f"   Keys: {list(adata.layers.keys())}")

    # uns
    print(f"\n7. Unstructured (adata.uns):")
    print(f"   Keys: {list(adata.uns.keys())}")

    # Save checkpoint
    output = Path("data/processed/rna/pbmc3k_inspected.h5ad")
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
    print(f"\nSaved to: {output}")

    print("\n" + "=" * 50)
    print("STEP 1 COMPLETE")
    print("=" * 50)

if __name__ == "__main__":
    main()
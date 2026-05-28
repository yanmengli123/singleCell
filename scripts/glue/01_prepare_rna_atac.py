#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 1: Prepare RNA and ATAC for GLUE integration
Uses scGLUE for graph-based multi-omics integration
"""
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def main():
    print("=" * 50)
    print("STEP 1: PREPARE GLUE INPUT")
    print("=" * 50)

    # Load GLUE tutorial data (Chen-2019)
    rna_file = Path("D:/singleCelldata/glue/Chen-2019-RNA.h5ad")
    atac_file = Path("D:/singleCelldata/glue/Chen-2019-ATAC.h5ad")
    output_dir = Path("data/processed/glue")

    print("\n[1] Loading GLUE tutorial data...")
    adata_rna = sc.read_h5ad(rna_file)
    adata_atac = sc.read_h5ad(atac_file)

    print(f"RNA: {adata_rna.n_obs} cells x {adata_rna.n_vars} genes")
    print(f"ATAC: {adata_atac.n_obs} cells x {adata_atac.n_vars} peaks")

    if 'cell_type' in adata_rna.obs.columns:
        print(f"\nRNA cell types: {adata_rna.obs['cell_type'].nunique()}")
        print(adata_rna.obs['cell_type'].value_counts())

    if 'cell_type' in adata_atac.obs.columns:
        print(f"\nATAC cell types: {adata_atac.obs['cell_type'].nunique()}")
        print(adata_atac.obs['cell_type'].value_counts())

    # Save for GLUE pipeline
    print("\n[2] Saving GLUE input data...")
    output_dir.mkdir(parents=True, exist_ok=True)
    adata_rna.write_h5ad(output_dir / "rna_prepared.h5ad")
    adata_atac.write_h5ad(output_dir / "atac_prepared.h5ad")

    print("\n" + "=" * 50)
    print("STEP 1 COMPLETE")
    print("=" * 50)

if __name__ == "__main__":
    main()
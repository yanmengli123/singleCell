#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 1: Prepare MultiVI input from Multiome data
Uses scvi-tools for VAE-based integration
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
    print("STEP 1: PREPARE MULTIVI INPUT")
    print("=" * 50)

    # Load processed multiome RNA and ATAC
    rna_file = Path("data/processed/multiome/rna_processed.h5ad")
    output_file = Path("data/processed/multivi/multivi_latent.h5ad")
    output_dir = Path("results/figures/multivi")

    print("\n[1] Loading processed multiome data...")
    adata_rna = sc.read_h5ad(rna_file)
    print(f"RNA: {adata_rna.n_obs} cells x {adata_rna.n_vars} genes")

    # For MultiVI, we need both RNA and ATAC
    # Since we have multiome data (paired), we can use both
    # MultiVI will create a joint latent space

    # Load ATAC data
    atac_file = Path("data/raw/multiome/atac.h5ad")
    adata_atac = sc.read_h5ad(atac_file)

    # Ensure cells match between RNA and ATAC
    common_cells = adata_rna.obs.index.intersection(adata_atac.obs.index)
    print(f"\n[2] Matching cells between RNA and ATAC...")
    print(f"  Common cells: {len(common_cells)}")

    adata_rna = adata_rna[common_cells].copy()
    adata_atac = adata_atac[common_cells].copy()

    # Save preprocessed data
    print("\n[3] Saving preprocessed data for MultiVI...")
    output_dir.mkdir(parents=True, exist_ok=True)
    adata_rna.write_h5ad("data/processed/multivi/rna_for_multivi.h5ad")
    adata_atac.write_h5ad("data/processed/multivi/atac_for_multivi.h5ad")

    print("\n" + "=" * 50)
    print("STEP 1 COMPLETE")
    print("=" * 50)

if __name__ == "__main__":
    main()
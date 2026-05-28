#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 1: Prepare Multiome (RNA + ATAC) data from 10x Genomics
"""
import scanpy as sc
import anndata as ad
import muon as mu
from muon import atac as atac
import numpy as np
import pandas as pd
from scipy import sparse
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


def read_10x_multiome_h5(filename):
    """Read 10x Multiome h5 file (RNA + ATAC)."""
    with h5py.File(filename, 'r') as f:
        print(f"Keys in h5: {list(f.keys())}")

        if 'matrix' in f:
            matrix = f['matrix']
            shape = matrix['shape'][:].tolist()  # [n_genes, n_cells]

            data = matrix['data'][:]
            indices = matrix['indices'][:]
            indptr = matrix['indptr'][:]

            # CSC matrix (n_genes, n_cells)
            csc = sparse.csc_matrix((data, indices, indptr), shape=tuple(shape))
            X = csc.T.tocsr()  # (n_cells, n_genes)

            barcodes = matrix['barcodes'][:].astype(str)
            features = matrix['features']

            # Get feature info
            feature_ids = features['id'][:].astype(str) if 'id' in features else None
            feature_names = features['name'][:].astype(str) if 'name' in features else None
            feature_types = features['feature_type'][:].astype(str) if 'feature_type' in features else None

            # Create AnnData
            adata = ad.AnnData(
                X=X,
                obs=pd.DataFrame(index=barcodes),
                var=pd.DataFrame({
                    'gene_ids': feature_ids,
                    'feature_types': feature_types
                }, index=feature_names)
            )

    return adata


def main():
    input_file = Path("D:/singleCelldata/multiome/filtered_feature_bc_matrix.h5")
    output_rna = Path("data/raw/multiome/rna.h5ad")
    output_atac = Path("data/raw/multiome/atac.h5ad")

    print("=" * 50)
    print("STEP 1: PREPARE MULTIOME DATA")
    print("=" * 50)

    import h5py

    # Read the multiome h5
    print("\n[1] Reading multiome h5 file...")
    adata = read_10x_multiome_h5(input_file)
    print(f"Combined data: {adata.n_obs} cells x {adata.n_vars} features")

    if 'feature_types' in adata.var.columns:
        print("\nFeature types:")
        print(adata.var['feature_types'].value_counts())

    # Split into RNA and ATAC
    print("\n[2] Splitting into RNA and ATAC...")
    rna_mask = adata.var['feature_types'] == 'Gene Expression'
    atac_mask = adata.var['feature_types'] == 'Peaks'

    print(f"  RNA features: {rna_mask.sum()}")
    print(f"  ATAC features: {atac_mask.sum()}")

    adata_rna = adata[:, rna_mask].copy()
    adata_atac = adata[:, atac_mask].copy()

    print(f"\nRNA: {adata_rna.n_obs} cells x {adata_rna.n_vars} genes")
    print(f"ATAC: {adata_atac.n_obs} cells x {adata_atac.n_vars} peaks")

    # Save
    print("\n[3] Saving RNA and ATAC separately...")
    output_rna.parent.mkdir(parents=True, exist_ok=True)
    adata_rna.write_h5ad(output_rna)
    print(f"  Saved RNA to: {output_rna}")

    adata_atac.write_h5ad(output_atac)
    print(f"  Saved ATAC to: {output_atac}")

    # Create MuData (combined multi-omics object)
    print("\n[4] Creating MuData object...")
    mdata = mu.MuData({
        'rna': adata_rna,
        'atac': adata_atac
    })
    print(mdata)

    print("\n" + "=" * 50)
    print("STEP 1 COMPLETE")
    print("=" * 50)


if __name__ == "__main__":
    import h5py
    main()
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 1: Inspect scATAC data structure
Reads PBMC 10k scATAC peak matrix from 10x Genomics
"""
import scanpy as sc
import anndata as ad
import numpy as np
import h5py
import pandas as pd
from scipy import sparse
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')


def read_10x_atac_h5(filename):
    """Read 10x scATAC h5 file (matrix format)."""
    with h5py.File(filename, 'r') as f:
        matrix = f['matrix']
        shape = matrix['shape'][:].tolist()  # [n_peaks, n_cells]
        data = matrix['data'][:]
        indices = matrix['indices'][:]
        indptr = matrix['indptr'][:]

        # Create CSC matrix (n_peaks, n_cells) then transpose
        csc = sparse.csc_matrix((data, indices, indptr), shape=tuple(shape))
        X = csc.T.tocsr()  # (n_cells, n_peaks)

        # Read barcodes
        barcodes = matrix['barcodes'][:].astype(str)

        # Read features (peaks)
        features = matrix['features']
        peak_names = features['name'][:].astype(str) if 'name' in features else [f'peak_{i}' for i in range(shape[0])]

        # Create AnnData
        adata = ad.AnnData(
            X=X,
            obs=pd.DataFrame(index=barcodes),
            var=pd.DataFrame(index=peak_names)
        )

    return adata


def main():
    input_file = Path("D:/singleCelldata/atac/filtered_peak_bc_matrix.h5")
    output_file = Path("data/raw/atac/pbmc10k_atac_raw.h5ad")
    output_dir = Path("results/figures/atac")

    print("=" * 50)
    print("STEP 1: INSPECT ATAC DATA")
    print("=" * 50)

    # Load 10x scATAC peak matrix
    print("\n[1] Loading scATAC peak matrix...")
    adata = read_10x_atac_h5(input_file)
    adata.var_names_make_unique()

    print(f"\nRaw data: {adata.n_obs} cells x {adata.n_vars} peaks")
    print(f"Matrix type: {type(adata.X)}")
    print(f"Matrix density: {adata.X.nnz / (adata.n_obs * adata.n_vars):.4%}")

    # Basic inspection
    print("\n[2] Inspecting data structure...")
    print(f"obs columns: {list(adata.obs.columns)}")
    print(f"var columns: {list(adata.var.columns)}")

    # Add basic QC metrics
    print("\n[3] Calculating basic QC...")
    adata.obs['n_peaks'] = np.array(adata.X.sum(axis=1)).flatten()
    adata.obs['total_counts'] = adata.obs['n_peaks'].values

    print(f"Median peaks per cell: {adata.obs['n_peaks'].median():.0f}")
    print(f"Mean peaks per cell: {adata.obs['n_peaks'].mean():.0f}")

    # Add cell metadata from singlecell.csv
    print("\n[4] Adding cell metadata...")
    sc_meta = pd.read_csv("D:/singleCelldata/atac/singlecell.csv", index_col=0)
    print(f"Metadata cells: {len(sc_meta)}")
    print(f"Metadata columns: {list(sc_meta.columns)[:10]}...")

    # Merge metadata
    common_cells = adata.obs.index.intersection(sc_meta.index)
    if len(common_cells) > 0:
        # Add columns that don't already exist
        for col in sc_meta.columns:
            if col not in adata.obs.columns:
                adata.obs[col] = sc_meta.loc[adata.obs.index, col]
        print(f"Matched cells: {len(common_cells)}")
    else:
        print("WARNING: No matching cells found")

    # Save raw ATAC data
    print("\n[5] Saving raw ATAC data...")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_file)
    print(f"  Saved to: {output_file}")

    # Quick QC plot
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    axes[0].hist(adata.obs['n_peaks'].values, bins=50, color='steelblue', edgecolor='white')
    axes[0].set_xlabel('Peaks Detected')
    axes[0].set_ylabel('Number of Cells')
    axes[0].set_title('Peaks per Cell')
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)

    if 'passed_filters' in adata.obs.columns:
        axes[1].hist(adata.obs['passed_filters'].values, bins=50, color='lightcoral', edgecolor='white')
        axes[1].set_xlabel('Total Fragments')
        axes[1].set_ylabel('Number of Cells')
        axes[1].set_title('Fragments per Cell')
        axes[1].spines['top'].set_visible(False)
        axes[1].spines['right'].set_visible(False)

    if 'duplicat_rate' in adata.obs.columns:
        dup_rate = adata.obs['duplicat_rate'].values
        axes[2].hist(dup_rate, bins=50, color='lightgreen', edgecolor='white')
        axes[2].set_xlabel('Duplicate Rate')
        axes[2].set_ylabel('Number of Cells')
        axes[2].set_title('Duplication Rate')
        axes[2].spines['top'].set_visible(False)
        axes[2].spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_dir / 'atac_qc_overview.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_dir / 'atac_qc_overview.png'}")

    # Summary statistics
    print("\n" + "=" * 50)
    print("SUMMARY:")
    print("=" * 50)
    print(f"Cells: {adata.n_obs}")
    print(f"Peaks: {adata.n_vars}")
    print(f"Median peaks/cell: {adata.obs['n_peaks'].median():.0f}")
    print(f"Total UMIs (median): {adata.obs['total_counts'].median():.0f}")

    print("\n" + "=" * 50)
    print("STEP 1 COMPLETE")
    print("=" * 50)

if __name__ == "__main__":
    main()
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ATAC Step 1: Inspect ATAC data structure.
This script is a placeholder - requires actual ATAC data.
"""
from pathlib import Path


def main():
    print("=" * 50)
    print("STEP 1: INSPECT ATAC DATA (PLACEHOLDER)")
    print("=" * 50)

    print("""
    ATAC data typically comes as:
    - cell x peak matrix (binary accessibility)
    - Fragment file (fragments.tsv.gz from 10x)
    - Peak annotation file

    Expected structure for muon/Signac processing:
    - adata_atac.X: sparse binary matrix
    - adata_atac.var: peak coordinates and annotations
    - adata_atac.obs: cell metadata

    To get test ATAC data:
    1. Download from 10x Genomics (PBMC scATAC)
    2. Or use: scanpy.datasets.pbmc3k_processed() for combined

    Placeholder script - implement when ATAC data is available.
    """)


if __name__ == "__main__":
    main()
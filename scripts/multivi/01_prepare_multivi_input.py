#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
MultiVI Step 1: Prepare MuData input for MultiVI.
"""
from pathlib import Path


def main():
    print("=" * 50)
    print("STEP 1: PREPARE MULTIVI INPUT (PLACEHOLDER)")
    print("=" * 50)

    print("""
    MultiVI uses scvi-tools for VAE-based integration.

    Input: MuData (.h5mu) object with:
    - mdata['rna']: RNA AnnData
    - mdata['atac']: ATAC AnnData

    Setup steps:
    1. Prepare RNA AnnData (normalized, log-transformed)
    2. Prepare ATAC AnnData (TF-IDF, binarized)
    3. Create MuData object
    4. Setup MultiVI model

    MultiVI handles:
    - Missing modalities (unpaired cells)
    - Batch effects
    - Joint latent space
    """)


if __name__ == "__main__":
    main()
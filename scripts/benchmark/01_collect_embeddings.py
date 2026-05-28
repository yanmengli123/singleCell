#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Benchmark Step 1: Collect embeddings from all methods.
"""
from pathlib import Path


def main():
    print("=" * 50)
    print("STEP 1: COLLECT EMBEDDINGS (PLACEHOLDER)")
    print("=" * 50)

    print("""
    Embedding collection for benchmark:

    Methods to compare:
    1. Scanpy baseline (RNA-only PCA)
    2. WNN (from multiome pipeline)
    3. MultiVI (from multivi pipeline)
    4. GLUE (from glue pipeline)

    For each method, collect:
    - adata.obsm['X_embedding'] or similar
    - adata.obs['cell_type']
    - adata.obs['batch'] (if available)
    - adata.obs['modality'] (RNA/ATAC)

    Save as:
    - data/processed/benchmark/{method}_embedding.h5ad
    """)


if __name__ == "__main__":
    main()
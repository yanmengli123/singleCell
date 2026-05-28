#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Benchmark Step 2: Compute metrics for each method.
"""
from pathlib import Path


def main():
    print("=" * 50)
    print("STEP 2: COMPUTE METRICS (PLACEHOLDER)")
    print("=" * 50)

    print("""
    Metrics to compute:

    Clustering quality:
    - ARI (Adjusted Rand Index)
    - NMI (Normalized Mutual Information)
    - silhouette score per cell type

    Batch correction:
    - silhouette score per batch (lower is better)
    - silhouette score per modality

    Modality mixing:
    - Modality mixing score (do modalities align?)
    - k-NN purity

    Performance:
    - Runtime (seconds)
    - GPU memory (GB)

    Comparison table structure:
    method | dataset | ARI | NMI | silhouette_celltype | ...
    """)


if __name__ == "__main__":
    main()
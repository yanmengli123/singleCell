#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ATAC Step 3: Gene activity matrix computation.
This script is a placeholder - requires actual ATAC data.
"""
from pathlib import Path


def main():
    print("=" * 50)
    print("STEP 3: GENE ACTIVITY (PLACEHOLDER)")
    print("=" * 50)

    print("""
    Gene activity matrices bridge ATAC and RNA:
    - Count peaks near genes
    - Weight by distance to TSS
    - Create pseudo-RNA expression matrix

    Uses:
    - sc.tl.score_gene_sets (Scanpy)
    - pycistopic topic scores
    - Signac gene activity matrix

    Output can be used for:
    - Cross-modality clustering
    - Multiome integration
    - Peak-to-gene linking

    Placeholder script - implement when ATAC data is available.
    """)


if __name__ == "__main__":
    main()
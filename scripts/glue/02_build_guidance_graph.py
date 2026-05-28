#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GLUE Step 2: Build guidance graph for RNA-ATAC integration.
"""
from pathlib import Path


def main():
    print("=" * 50)
    print("STEP 2: BUILD GUIDANCE GRAPH (PLACEHOLDER)")
    print("=" * 50)

    print("""
    Guidance graph connects RNA features to ATAC features:

    Node types:
    - Gene nodes (RNA features)
    - Peak nodes (ATAC features)
    - Optional: TF nodes, Motif nodes

    Edge types:
    - Gene-Peak: based on genomic proximity (distance to TSS)
    - TF-Motif: transcription factor binding motifs
    - Correlation: co-accessibility

    Graph building steps:
    1. Select HVG for RNA
    2. Select highly variable peaks for ATAC
    3. Link peaks near gene TSS
    4. Filter by correlation or distance threshold
    5. Export as GraphML for scGLUE

    Guidance graph provides supervision signal
    that helps align RNA and ATAC embeddings.
    """)


if __name__ == "__main__":
    main()
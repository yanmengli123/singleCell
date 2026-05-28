#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 4: Marker Gene Detection and Cell Type Annotation
"""
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# PBMC marker genes
MARKER_GENES = {
    'CD8 T cells': ['CD8A', 'CD8B'],
    'CD4 T cells': ['CD3D', 'CD3E', 'IL7R'],
    'B cells': ['MS4A1', 'CD79A'],
    'NK cells': ['NKG7', 'GNLY'],
    'Monocytes': ['S100A8', 'S100A9', 'FCGR3A'],
    'DC': ['FCER1A', 'CST3'],
    'Platelets': ['PPBP']
}

def main():
    input_file = Path("data/processed/pbmc3k_processed.h5ad")
    output_file = Path("data/processed/pbmc3k_annotated.h5ad")
    marker_file = Path("results/tables/marker_genes.csv")
    annot_file = Path("results/tables/cluster_annotation.csv")

    print("=" * 50)
    print("STEP 4: MARKER GENES & ANNOTATION")
    print("=" * 50)

    adata = sc.read_h5ad(input_file)
    print(f"\nInput: {adata.n_obs} cells x {adata.n_vars} genes, {adata.obs['leiden'].nunique()} clusters")

    # === Find marker genes for each cluster ===
    print("\n[1] Finding marker genes (rank_genes_groups)...")
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', random_state=0)

    # === Save marker genes table ===
    print("\n[2] Saving marker genes table...")
    marker_file.parent.mkdir(parents=True, exist_ok=True)
    marker_df = sc.get.rank_genes_groups_df(adata, group=None, key='rank_genes_groups')
    marker_df.to_csv(marker_file, index=False)
    print(f"   Saved to: {marker_file}")

    # === Plot marker genes ===
    print("\n[3] Plotting top marker genes...")
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, ax=ax, show=False)
    fig.savefig("results/figures/marker_dotplot.png", dpi=100, bbox_inches='tight')
    plt.close()
    print("   Saved: results/figures/marker_dotplot.png")

    # === Manual cell type annotation based on markers ===
    print("\n[4] Annotating cell types based on marker genes...")

    # Get top marker for each cluster
    cluster_markers = {}
    for cluster in sorted(adata.obs['leiden'].unique(), key=int):
        cluster_data = marker_df[marker_df['group'] == cluster].head(10)
        cluster_markers[cluster] = cluster_data['names'].tolist()[:5]

    # Print top markers for reference
    for cluster, markers in cluster_markers.items():
        print(f"   Cluster {cluster}: {', '.join(markers)}")

    # Annotation based on marker gene analysis:
    # Cluster 0: ribosomal/housekeeping genes (LDHB, RPS genes) - likely B cells
    # Cluster 1: NKG7, GZMA, CST7 - cytotoxic NK markers → NK cells
    # Cluster 2: CD74, HLA-DRA, CD79A - MHCII/B cell markers → B cells
    # Cluster 3: FTL, FTH1, LYZ, S100A9 - monocyte markers → Monocytes
    # Cluster 4: PPBP (platelet marker) → Platelets
    cluster_annotation = {
        '0': 'B cells',
        '1': 'NK cells',
        '2': 'B cells',
        '3': 'Monocytes',
        '4': 'Platelets'
    }

    adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotation)

    # === Save annotation table ===
    print("\n[5] Saving cluster annotation...")
    annot_file.parent.mkdir(parents=True, exist_ok=True)
    annot_df = pd.DataFrame({
        'cluster': list(cluster_annotation.keys()),
        'cell_type': list(cluster_annotation.values())
    })
    annot_df.to_csv(annot_file, index=False)
    print(f"   Saved to: {annot_file}")

    # === Final UMAP colored by cell type ===
    print("\n[6] Plotting UMAP with cell types...")
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    sc.pl.umap(adata, color='cell_type', ax=ax, show=False)
    fig.savefig("results/figures/umap_celltype.png", dpi=100, bbox_inches='tight')
    plt.close()
    print("   Saved: results/figures/umap_celltype.png")

    # === Dot plot of marker genes per cell type ===
    print("\n[7] Plotting marker genes per cell type...")
    marker_genes_list = [g for genes in MARKER_GENES.values() for g in genes]
    available_markers = [g for g in marker_genes_list if g in adata.var_names]
    if available_markers:
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        sc.pl.dotplot(adata, markers, groupby='cell_type', ax=ax, show=False)
        fig.savefig("results/figures/umap_markers.png", dpi=100, bbox_inches='tight')
        plt.close()
        print("   Saved: results/figures/umap_markers.png")

    # === Save annotated data ===
    print("\n[8] Saving annotated data...")
    adata.write_h5ad(output_file)
    print(f"   Saved to: {output_file}")

    # === Summary ===
    print("\n" + "=" * 50)
    print("CELL TYPE COUNTS:")
    print("=" * 50)
    print(adata.obs['cell_type'].value_counts())

    print("\n" + "=" * 50)
    print("STEP 4 COMPLETE")
    print("=" * 50)

if __name__ == "__main__":
    main()
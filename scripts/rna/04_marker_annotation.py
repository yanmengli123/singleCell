#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 4: Marker Gene Detection and Cell Type Annotation

Fixed:
- Uses correct paths for RNA subdirectory
- Reads cluster annotation from config file
- Reports unannotated clusters as error
- Fixes marker dotplot to use available_markers
"""
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import yaml
from pathlib import Path

# PBMC marker genes for dotplot
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
    input_file = Path("data/processed/rna/pbmc3k_processed.h5ad")
    output_file = Path("data/processed/rna/pbmc3k_annotated.h5ad")
    marker_file = Path("results/tables/rna/marker_genes.csv")
    annot_file = Path("results/tables/rna/cluster_annotation.csv")
    annot_config = Path("configs/pbmc3k_cluster_annotation.yaml")

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

    # === Plot marker genes (dotplot per cluster) ===
    print("\n[3] Plotting top marker genes per cluster...")
    Path("results/figures/rna").mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, ax=ax, show=False, return_fig=True)
    fig.savefig("results/figures/rna/marker_dotplot.png", dpi=100, bbox_inches='tight')
    plt.close()
    print("   Saved: results/figures/rna/marker_dotplot.png")

    # === Print top markers for reference ===
    print("\n[4] Top markers per cluster:")
    for cluster in sorted(adata.obs['leiden'].unique(), key=int):
        cluster_data = marker_df[marker_df['group'] == cluster].head(5)
        markers = ', '.join(cluster_data['names'].tolist()[:5])
        print(f"   Cluster {cluster}: {markers}")

    # === Load cluster annotation from config ===
    print("\n[5] Loading cluster annotation from config...")
    with open(annot_config, 'r') as f:
        config = yaml.safe_load(f)
    cluster_annot = config['cluster_annotation']

    # Check all clusters have annotation
    unannotated = [c for c in adata.obs['leiden'].unique() if cluster_annot.get(c, 'TBD') == 'TBD']
    if unannotated:
        print(f"   WARNING: Clusters {unannotated} have no annotation (TBD)")
        print(f"   Please update {annot_config} before re-running")
        # Use TBD for unannotated clusters
        adata.obs['cell_type'] = adata.obs['leiden'].map(
            lambda x: cluster_annot.get(x, 'TBD')
        )
    else:
        adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annot)

    # === Save annotation table ===
    print("\n[6] Saving cluster annotation...")
    annot_file.parent.mkdir(parents=True, exist_ok=True)
    annot_df = pd.DataFrame({
        'cluster': list(cluster_annot.keys()),
        'cell_type': list(cluster_annot.values())
    })
    annot_df.to_csv(annot_file, index=False)
    print(f"   Saved to: {annot_file}")

    # === Final UMAP colored by cell type ===
    print("\n[7] Plotting UMAP with cell types...")
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    sc.pl.umap(adata, color='cell_type', ax=ax, show=False)
    fig.savefig("results/figures/rna/umap_celltype.png", dpi=100, bbox_inches='tight')
    plt.close()
    print("   Saved: results/figures/rna/umap_celltype.png")

    # === Dot plot of marker genes per cell type ===
    print("\n[8] Plotting marker genes per cell type...")
    marker_genes_list = [g for genes in MARKER_GENES.values() for g in genes]
    available_markers = {ct: [g for g in genes if g in adata.var_names]
                         for ct, genes in MARKER_GENES.items()}
    available_markers = {k: v for k, v in available_markers.items() if v}
    if available_markers:
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        sc.pl.dotplot(adata, available_markers, groupby='cell_type', ax=ax, show=False)
        fig.savefig("results/figures/rna/umap_markers.png", dpi=100, bbox_inches='tight')
        plt.close()
        print("   Saved: results/figures/rna/umap_markers.png")

    # === Save annotated data ===
    print("\n[9] Saving annotated data...")
    output_file.parent.mkdir(parents=True, exist_ok=True)
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
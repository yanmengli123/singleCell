# singleCell

PBMC scRNA-seq standard analysis project using Scanpy.

## Project Structure

```
singleCell/
  data/
    raw/           # Raw h5ad files
      pbmc3k_raw.h5ad
    processed/     # Processed h5ad files
      pbmc3k_inspected.h5ad
      pbmc3k_qc.h5ad
      pbmc3k_processed.h5ad
      pbmc3k_annotated.h5ad
  scripts/         # Analysis scripts
    00_check_environment.py
    01_inspect_anndata.py
    02_qc_filter.py
    03_preprocess_cluster.py
    04_marker_annotation.py
  notebooks/        # Jupyter notebooks
    01_pbmc3k_rna_scanpy.ipynb
  results/
    figures/        # Output plots
    tables/         # Output CSVs
  notes/           # Project notes
  configs/         # Configuration files
  logs/            # Run logs
```

## Quick Start

```bash
# Check environment
D:/soft/Python310/python.exe scripts/00_check_environment.py

# Run full pipeline
D:/soft/Python310/python.exe scripts/01_inspect_anndata.py
D:/soft/Python310/python.exe scripts/02_qc_filter.py
D:/soft/Python310/python.exe scripts/03_preprocess_cluster.py
D:/soft/Python310/python.exe scripts/04_marker_annotation.py
```

## Pipeline Summary

| Step | Input | Output |
|------|-------|--------|
| 1. Inspect | raw h5ad | inspected h5ad |
| 2. QC | inspected | qc filtered |
| 3. Cluster | qc | processed + PCA/UMAP/Leiden |
| 4. Annotate | processed | annotated + cell types |

## Results

- Cells: 2700 → 1705 (after QC)
- Genes: 32738 → 12360 (after filtering)
- Clusters: 5
- Cell types: B cells, NK cells, Monocytes, Platelets
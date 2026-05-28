# singleCell

Single-cell multiomics analysis project: RNA, scATAC, 10x Multiome, MultiVI, GLUE, and benchmark workflows.

## Current Status

**Operational**: PBMC3k RNA-only pipeline
**Placeholder**: scATAC, Multiome, MultiVI, GLUE, Benchmark (require additional data)

## Python Interpreter

```
D:/soft/Python310/python.exe
```

## Quick Start - RNA Pipeline

```bash
# Check environment
D:/soft/Python310/python.exe scripts/rna/00_check_environment.py

# Run pipeline steps in order
D:/soft/Python310/python.exe scripts/rna/01_inspect_anndata.py
D:/soft/Python310/python.exe scripts/rna/02_qc_filter.py
D:/soft/Python310/python.exe scripts/rna/03_preprocess_cluster.py
D:/soft/Python310/python.exe scripts/rna/04_marker_annotation.py
```

## Pipeline Summary

| Step | Input | Output |
|------|-------|--------|
| 1. Inspect | raw h5ad | inspected h5ad |
| 2. QC | raw | qc filtered (cells+genes) |
| 3. Cluster | qc | processed + PCA/UMAP/Leiden |
| 4. Annotate | processed | annotated + cell types |

## Results (After Fix)

- Cells: 2700 -> 2638 (after QC with correct n_genes_by_counts)
- Genes: 32738 -> 13656 (after filtering)
- Clusters: 6 (at resolution 0.5)
- Cell types: CD4 T cells, B cells, Monocytes, NK cells, Non-classical Monocytes, Dendritic cells
- Raw counts preserved in `adata.layers['counts']`

## Project Structure

```
singleCell/
  configs/              # YAML configs for all parameters
  scripts/
    common/            # Shared utilities
    rna/                # RNA pipeline (operational)
    atac/               # ATAC pipeline (placeholder)
    multiome/           # Multiome pipeline (placeholder)
    multivi/            # MultiVI pipeline (placeholder)
    glue/               # GLUE pipeline (placeholder)
    benchmark/          # Benchmark scripts (placeholder)
  data/
    raw/rna/            # Raw data (read-only, tracked)
    processed/rna/       # Intermediate files (regenerate via scripts)
  results/
    figures/rna/        # Output plots
    tables/rna/         # Output CSVs
  docs/                 # Project documentation
  tests/                # Data integrity tests
```

## Key Design

- **QC metrics**: n_genes_by_counts = number of genes with >0 counts (not sum of counts)
- **Paths**: All paths from configs/*.yaml, not hardcoded
- **Layers**: Raw counts preserved in `adata.layers['counts']`
- **Random state**: All stochastic ops use `random_state=0`
- **Annotation**: Cluster annotations read from configs/pbmc3k_cluster_annotation.yaml

## Large File Policy

- `data/processed/*.h5ad` >100MB: git-ignored, regenerate with scripts
- `results/figures/` and `results/tables/`: git-ignored, regenerate from scripts
- Original data in `data/raw/rna/` tracked in git
# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Single-cell multiomics analysis project. **RNA-only pipeline (PBMC3k) is operational.** Other modalities (scATAC, Multiome, MultiVI, GLUE, Benchmark) are placeholder scripts.

## Python Interpreter

```
D:/soft/Python310/python.exe
```

## Running the RNA Pipeline

```bash
# Full pipeline (run in order)
D:/soft/Python310/python.exe scripts/rna/00_check_environment.py
D:/soft/Python310/python.exe scripts/rna/01_inspect_anndata.py
D:/soft/Python310/python.exe scripts/rna/02_qc_filter.py
D:/soft/Python310/python.exe scripts/rna/03_preprocess_cluster.py
D:/soft/Python310/python.exe scripts/rna/04_marker_annotation.py
```

## Running Tests

```bash
D:/soft/Python310/python.exe tests/test_data_integrity.py
```

Tests validate: file existence, config paths match scripts, QC correctness (n_genes_by_counts != total_counts), cell_type has no NA, all clusters annotated.

## Architecture

### Pipeline Structure
- **scripts/rna/**: 5-step RNA pipeline (operational)
- **scripts/common/**: Shared utilities (qc_utils.py, plotting_utils.py, etc.)
- **scripts/atac/, multiome/, multivi/, glue/, benchmark/**: Placeholder scripts
- **configs/**: YAML configs with all parameters (not hardcoded)
- **tests/**: Data integrity tests

### Data Flow
- `data/raw/rna/pbmc3k_raw.h5ad` (read-only)
- `data/processed/rna/` intermediate files (>100MB git-ignored)
- `results/figures/rna/` and `results/tables/rna/` (git-ignored)

### Key Design Patterns

**QC Metrics** (`scripts/common/qc_utils.py`):
- `n_genes_by_counts` = number of genes with >0 counts per cell (NOT sum of counts)
- `total_counts` = sum of all counts per cell
- `pct_counts_mt` = mitochondrial percentage

**Layers**: Raw counts preserved in `adata.layers['counts']` before normalization.

**Cluster Annotation**: Read from `configs/pbmc3k_cluster_annotation.yaml`, not hardcoded in scripts.

**Random State**: All stochastic operations use `random_state=0` for reproducibility.

## Important Notes

- `data/processed/*.h5ad` >100MB are git-ignored; regenerate from scripts
- Figures/tables in `results/` are git-ignored; regenerate from scripts
- Original data `data/raw/rna/pbmc3k_raw.h5ad` is tracked in git
- Console output: use UTF-8 compatible strings (avoid emoji/special chars)
- Paths: read from configs, do not hardcode in scripts
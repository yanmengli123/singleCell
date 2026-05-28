# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Single-cell multiomics analysis project. **All major pipelines are now operational** (RNA, scATAC, Multiome, GLUE, Benchmark).

## Python Interpreter

```
D:/soft/Python310/python.exe
```

## Running Pipelines

### RNA Pipeline (PBMC3k)
```bash
D:/soft/Python310/python.exe scripts/rna/00_check_environment.py
D:/soft/Python310/python.exe scripts/rna/01_inspect_anndata.py
D:/soft/Python310/python.exe scripts/rna/02_qc_filter.py
D:/soft/Python310/python.exe scripts/rna/03_preprocess_cluster.py
D:/soft/Python310/python.exe scripts/rna/04_marker_annotation.py
D:/soft/Python310/python.exe scripts/rna/05_visualization.py
```

### scATAC Pipeline (PBMC 10k)
```bash
D:/soft/Python310/python.exe scripts/atac/01_inspect_atac.py
D:/soft/Python310/python.exe scripts/atac/02_qc_lsi.py
```

### Multiome Pipeline
```bash
D:/soft/Python310/python.exe scripts/multiome/01_prepare_multiome.py
D:/soft/Python310/python.exe scripts/multiome/02_rna_only.py
```

### GLUE Pipeline
```bash
D:/soft/Python310/python.exe scripts/glue/01_prepare_rna_atac.py
```

### Benchmark
```bash
D:/soft/Python310/python.exe scripts/benchmark/01_collect_embeddings.py
```

## Running Tests

```bash
D:/soft/Python310/python.exe tests/test_data_integrity.py
```

Tests validate: file existence, config paths match scripts, QC correctness (n_genes_by_counts != total_counts), cell_type has no NA, all clusters annotated.

## Architecture

### Pipeline Structure
- **scripts/rna/**: 5-step RNA pipeline (operational)
- **scripts/atac/**: scATAC with TF-IDF + LSI (operational)
- **scripts/multiome/**: 10x Multiome RNA+ATAC (operational)
- **scripts/glue/**: Graph-based integration (operational)
- **scripts/benchmark/**: Method comparison (operational)
- **scripts/common/**: Shared utilities (qc_utils.py, plotting_utils.py, etc.)
- **configs/**: YAML configs with all parameters (not hardcoded)
- **tests/**: Data integrity tests

### External Data
Data is stored in `D:/singleCelldata/` (outside repo) to avoid large files in git:
- `D:/singleCelldata/atac/` - PBMC 10k scATAC (peak matrix + fragments)
- `D:/singleCelldata/multiome/` - PBMC 10k Multiome (RNA + ATAC)
- `D:/singleCelldata/multivi/` - MultiVI tutorial data (h5mu)
- `D:/singleCelldata/glue/` - GLUE Chen-2019 tutorial data

### Data Flow
- `data/raw/{modality}/` - Raw h5ad files (read-only)
- `data/processed/{modality}/` - Intermediate files (>100MB git-ignored)
- `results/figures/{modality}/` - Publication figures (300 DPI)
- `results/tables/{modality}/` - CSV results (git-ignored)

### Key Design Patterns

**QC Metrics** (`scripts/common/qc_utils.py`):
- `n_genes_by_counts` = number of genes with >0 counts per cell (NOT sum of counts)
- `total_counts` = sum of all counts per cell
- `pct_counts_mt` = mitochondrial percentage
- `n_peaks` = number of peaks detected per cell (ATAC)

**Layers**: Raw counts preserved in `adata.layers['counts']` before normalization.

**Cluster Annotation**: Read from `configs/pbmc3k_cluster_annotation.yaml`, not hardcoded in scripts.

**Random State**: All stochastic operations use `random_state=0` for reproducibility.

**Publication Figures**: All figures generated at 300 DPI using Arial font, frameless style.

## Important Notes

- `data/processed/*.h5ad` >100MB are git-ignored; regenerate from scripts
- Figures/tables in `results/` are git-ignored; regenerate from scripts
- Original data `data/raw/rna/pbmc3k_raw.h5ad` is tracked in git
- Console output: use UTF-8 compatible strings (avoid emoji/special chars)
- Paths: read from configs, do not hardcode in scripts
- scvi-tools and scglue cannot be installed on this Windows system (long path issue)
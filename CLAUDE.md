# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Single-cell multiomics analysis project spanning RNA-only, scATAC, 10x Multiome, MultiVI, GLUE, and benchmark workflows. Currently Phase 1: PBMC3k RNA-only pipeline is operational.

## Python Interpreter

Always use: `D:/soft/Python310/python.exe`

## Running the RNA Pipeline

```bash
# Check environment
D:/soft/Python310/python.exe scripts/rna/00_check_environment.py

# Run pipeline steps in order
D:/soft/Python310/python.exe scripts/rna/01_inspect_anndata.py
D:/soft/Python310/python.exe scripts/rna/02_qc_filter.py
D:/soft/Python310/python.exe scripts/rna/03_preprocess_cluster.py
D:/soft/Python310/python.exe scripts/rna/04_marker_annotation.py
```

## Architecture

### Pipeline Stages (8 phases total)
1. **RNA-only** (`scripts/rna/`) - Operational, uses PBMC3k data
2. **scATAC** (`scripts/atac/`) - Placeholder scripts, needs ATAC data
3. **Multiome + WNN** (`scripts/multiome/`) - Placeholder scripts, needs multiome data
4. **MultiVI** (`scripts/multivi/`) - VAE-based integration, placeholder
5. **GLUE** (`scripts/glue/`) - Graph-based integration, placeholder
6. **Benchmark** (`scripts/benchmark/`) - Method comparison, placeholder

### Data Flow
- Raw data: `data/raw/` (read-only, tracked in git)
- Processed: `data/processed/{modality}/` (regenerate via scripts, >100MB excluded from git)
- Results: `results/figures/` and `results/tables/` (regenerate from scripts)

### Configuration-Driven Design
All parameters are in YAML configs under `configs/`. Each pipeline stage has its own config:
- `01_pbmc3k_rna.yaml` - QC thresholds, HVG params, cell markers, output paths
- `05_glue.yaml` - Guidance graph, model architecture, training params

Modify configs to experiment; don't hardcode parameters in scripts.

### Common Utilities (`scripts/common/`)
- `io_utils.py` - File I/O helpers
- `qc_utils.py` - QC metric calculation and filtering
- `plotting_utils.py` - Figure saving with auto directory creation
- `logging_utils.py` - Pipeline step logging with timestamps
- `validation_utils.py` - AnnData shape and field validation

### Key Design Patterns
1. **Validation on load** - Scripts validate h5ad files exist and have required fields before processing
2. **Intermediate saves** - Each step saves its output for checkpointing
3. **Layer preservation** - Raw counts stored in `adata.layers['counts']` before normalization
4. **Random state** - All stochastic operations use `random_state=0` for reproducibility

## Important Notes

- Large h5ad files (>100MB) in `data/processed/` are git-ignored - regenerate with scripts
- Figures and tables are git-ignored - regenerate from scripts
- The original PBMC3k data is cached by scanpy at `~/.cache/scanpy/`
- Console output uses UTF-8; avoid emoji characters that cause GBK encoding errors on Windows
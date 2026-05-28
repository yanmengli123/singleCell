# Project Plan

## Overview
Single-cell multiomics learning, reproduction, and benchmark project starting from PBMC RNA-only and extending to scATAC, 10x Multiome, MultiVI, GLUE, and benchmark.

## Project Phases

### Phase 0: Environment Setup
- [x] Directory structure created
- [x] Common utilities implemented
- [x] Config files created
- [ ] Environment check script verified

### Phase 1: PBMC3k RNA-only (Week 1)
- [x] Raw data: 2700 cells x 32738 genes
- [x] QC filtering: 1705 cells x 12360 genes
- [x] PCA, Leiden, UMAP, cell type annotation
- [ ] Create report
- [ ] Tune parameters

### Phase 2: scATAC Basic (Week 2-3)
- [ ] Download test ATAC data
- [ ] Implement ATAC QC + LSI
- [ ] Gene activity matrix
- [ ] Create report

### Phase 3: 10x Multiome + WNN (Week 3-4)
- [ ] Prepare multiome inputs
- [ ] RNA-only pipeline
- [ ] ATAC-only pipeline
- [ ] WNN integration (R/Seurat)
- [ ] Create report

### Phase 4: MultiVI (Week 5)
- [ ] Prepare MuData input
- [ ] Train MultiVI
- [ ] Extract latent embedding
- [ ] Evaluate modality mixing

### Phase 5: GLUE (Week 6-7)
- [ ] Prepare RNA + ATAC inputs
- [ ] Build guidance graph
- [ ] Train GLUE
- [ ] Evaluate cell type transfer
- [ ] Ablation experiments

### Phase 6: Benchmark (Week 7-8)
- [ ] Collect embeddings from all methods
- [ ] Compute metrics (ARI, NMI, silhouette)
- [ ] Create summary tables
- [ ] Create report

### Phase 7: Paper Reproduction (Ongoing)
- [ ] Select paper from reading queue
- [ ] Reproduce one figure/table
- [ ] Document findings

## Directory Structure
```
singleCell/
├── data/
│   ├── raw/              # Original data (read-only)
│   ├── processed/        # Intermediate h5ad/h5mu files
│   └── external/          # References, annotations
├── scripts/
│   ├── common/           # Shared utilities
│   ├── rna/              # RNA pipeline
│   ├── atac/             # ATAC pipeline
│   ├── multiome/         # Multiome pipeline
│   ├── multivi/          # MultiVI pipeline
│   ├── glue/             # GLUE pipeline
│   └── benchmark/        # Benchmark scripts
├── configs/              # YAML config files
├── notebooks/            # Jupyter notebooks
├── results/              # Figures, tables, reports
├── docs/                 # Project documentation
├── logs/                 # Pipeline logs
└── tests/                # Data integrity tests
```

## Running the Pipeline

### RNA Pipeline
```bash
D:/soft/Python310/python.exe scripts/rna/00_check_environment.py
D:/soft/Python310/python.exe scripts/rna/01_inspect_anndata.py
D:/soft/Python310/python.exe scripts/rna/02_qc_filter.py
D:/soft/Python310/python.exe scripts/rna/03_preprocess_cluster.py
D:/soft/Python310/python.exe scripts/rna/04_marker_annotation.py
```

## Key Design Principles
1. Original data is read-only
2. Each step outputs intermediate h5ad/h5mu
3. All parameters in YAML configs
4. Logging for every step
5. Validation after each step
6. Modular, reusable code in scripts/common/
# Experiment Notes

## Week 1: PBMC3k RNA-only

### 2026-05-28
- Downloaded and loaded PBMC3k dataset from 10x Genomics
- Dataset: 2700 cells x 32738 genes
- Initial QC metrics calculated:
  - n_genes_by_counts
  - total_counts
  - pct_counts_mt (mitochondrial percentage)
- Filtering results:
  - Cells: 2700 -> 1705 (removed 995, 37%)
    - Most removed due to max_genes_per_cell=2500 threshold
    - 57 cells had high mitochondrial %
  - Genes: 32738 -> 12360 (removed 20378, 62%)
    - Filtered by min_cells=3
- HVG selected: 2000 genes
- PCA: 30 components
- Leiden clustering: 5 clusters at resolution=0.5
- Cell type annotation based on marker genes:
  - Cluster 0: B cells (ribosomal/housekeeping markers)
  - Cluster 1: NK cells (NKG7, GZMA)
  - Cluster 2: B cells (CD74, HLA-DRA, CD79A)
  - Cluster 3: Monocytes (FTL, FTH1, LYZ, S100A9)
  - Cluster 4: Platelets (PPBP)

### Issues Encountered
- Chinese console output encoding issue (GBK vs UTF-8)
  - Resolved by removing Chinese characters from print statements
- torch.cuda warning on import
  - pynvml deprecation warning, non-critical
- Sparse matrix warning during scale
  - Normal behavior for sparse data

### Next Steps
- [ ] Re-run with different Leiden resolution (0.3, 0.5, 0.8)
- [ ] Add CD4/CD8 T cell markers
- [ ] Generate QC violin plots
- [ ] Create week 1 report
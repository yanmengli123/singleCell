# Experiment Notes

## Week 1: PBMC3k RNA-only (Re-run after QC fix)

### 2026-05-28
- Downloaded and loaded PBMC3k dataset from 10x Genomics
- Dataset: 2700 cells x 32738 genes
- Initial QC metrics calculated (FIXED: n_genes_by_counts now correct):
  - n_genes_by_counts: number of genes with >0 counts (median 817)
  - total_counts: sum of all counts (median 2197)
  - pct_counts_mt (mitochondrial percentage)
- Filtering results (AFTER FIX):
  - Cells: 2700 -> 2638 (removed 62, 2.3%)
    - Only 5 cells removed due to max_genes_per_cell=2500 (was 942 before!)
    - 57 cells removed due to high mitochondrial %
  - Genes: 32738 -> 13656 (removed 19082, 58%)
  - Filtered by min_cells=3
- HVG selected: 2000 genes
- PCA: 30 components
- Leiden clustering: 6 clusters at resolution=0.5
- Cell type annotation based on marker genes:
  - Cluster 0: CD4 T cells (LDHB, RPS genes, CD3D)
  - Cluster 1: B cells (CD74, CD79A, MS4A1, HLA-DR)
  - Cluster 2: Monocytes (LYZ, S100A8, S100A9)
  - Cluster 3: NK cells (NKG7, GZMA, PRF1)
  - Cluster 4: Non-classical Monocytes (FCER1G, FCGR3A, LST1)
  - Cluster 5: Dendritic cells (FCER1A, HLA-DR)

### Issues Encountered
1. **CRITICAL QC BUG**: n_genes_by_counts was incorrectly calculated as total_counts
   - Fixed by changing `(adata.X.sum(axis=1))` to `(adata.X > 0).sum(axis=1)`
2. Chinese console output encoding issue (GBK vs UTF-8)
   - Resolved by removing Chinese characters from print statements
3. torch.cuda warning on import
   - pynvml deprecation warning, non-critical
4. Sparse matrix warning during scale
   - Normal behavior for sparse data

### Code Review Findings (Addressed)
1. n_genes_by_counts calculation - FIXED
2. Script paths - FIXED (now use data/raw/rna/ and data/processed/rna/)
3. Hardcoded cell annotation - FIXED (now reads from configs/pbmc3k_cluster_annotation.yaml)
4. PCA variance plot - FIXED (now uses adata.uns['pca']['variance_ratio'])
5. Marker dotplot variable - FIXED (now uses available_markers dict)
6. README paths - FIXED (now references scripts/rna/)

### Next Steps
- [x] Fix QC calculation (done - n_genes_by_counts now correct)
- [x] Fix script paths (done - now consistent with project structure)
- [x] Re-run pipeline (done - 2638 cells, 6 clusters)
- [ ] Generate QC violin plots
- [ ] Re-run with different Leiden resolution (0.3, 0.5, 0.8)
- [ ] Create week 1 report
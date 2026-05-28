#!/usr/bin/env python
# ==========================================
# 项目 1: PBMC scRNA-seq 标准流程
# Step 1: 数据加载与结构解剖
# ==========================================

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

# 设置 Scanpy 的显示参数
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# ==========================================
# 1. 下载并加载 PBMC 3k 数据集
# ==========================================
print("正在从 10x Genomics 下载 PBMC 3k 数据集...")
adata = sc.datasets.pbmc3k()
print("下载完成！")

# ==========================================
# 2. 解剖 AnnData 对象
# ==========================================
print("\n" + "="*30)
print("AnnData 对象概览")
print("="*30)
print(adata)

print("\n细胞元数据 (adata.obs) 前 5 行:")
print(adata.obs.head())

print("\n基因元数据 (adata.var) 前 5 行:")
print(adata.var.head())

print("\n表达矩阵 (adata.X) 信息:")
print(f"矩阵形状 (细胞数 x 基因数): {adata.shape}")
print(f"矩阵数据类型: {type(adata.X)}")

# ==========================================
# 3. 关键抽屉详解
# ==========================================
print("\n" + "="*30)
print("AnnData 核心抽屉详解")
print("="*30)
print(f"""
1. adata.X  - 表达矩阵 (Sparse matrix)
   形状: {adata.X.shape}
   数据类型: {type(adata.X)}

2. adata.obs - 细胞元数据 (DataFrame)
   列名: {list(adata.obs.columns)}
   行数 (细胞数): {adata.n_obs}

3. adata.var - 基因元数据 (DataFrame)
   列名: {list(adata.var.columns)}
   行数 (基因数): {adata.n_vars}

4. adata.obsm - 多维矩阵 (dict)
   包含的key: {list(adata.obsm.keys()) if hasattr(adata.obsm, 'keys') else '无'}
""")
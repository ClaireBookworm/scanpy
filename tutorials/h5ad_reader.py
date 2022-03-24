import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix

# print(ad.__version__)


sc.read_h5ad("data/pbmc3k_processed.h5ad")
adata = sc.read_h5ad("data/pbmc3k_processed.h5ad")

adata.X

adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
# print(adata.obs_names[:10])
# print (adata.var_names[:10])

# output
# Index(['Cell_0', 'Cell_1', 'Cell_2', 'Cell_3', 'Cell_4', 'Cell_5', 'Cell_6',
#        'Cell_7', 'Cell_8', 'Cell_9'],
#       dtype='object')
# Index(['Gene_0', 'Gene_1', 'Gene_2', 'Gene_3', 'Gene_4', 'Gene_5', 'Gene_6',
#        'Gene_7', 'Gene_8', 'Gene_9'],
#       dtype='object')

# rdata = ad.read('data/pbmc3k_processed.h5ad', backed='r')
# rdata.isbacked

print(adata[:5, ['Gene_1', 'Gene_3']])
# View of AnnData object with n_obs × n_vars = 5 × 2
#     obs: 'n_genes', 'percent_mito', 'n_counts', 'louvain'
#     var: 'n_cells'
#     uns: 'draw_graph', 'louvain', 'louvain_colors', 'neighbors', 'pca', 'rank_genes_groups'
#     obsm: 'X_pca', 'X_tsne', 'X_umap', 'X_draw_graph_fr'
#     varm: 'PCs'
#     obsp: 'distances', 'connectivities'
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
print(adata.obs_names[:10])

# output
# Index(['Cell_0', 'Cell_1', 'Cell_2', 'Cell_3', 'Cell_4', 'Cell_5', 'Cell_6',
#        'Cell_7', 'Cell_8', 'Cell_9'],
#       dtype='object')

# adata = ad.read('my_results.h5ad', backed='r')
# adata.isbacked
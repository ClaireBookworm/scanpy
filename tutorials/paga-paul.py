# paga paul15

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc

sc.settings.verbosity = 3
# sc.logging.print_versions()
results_file = "./write/paul15_results.h5ad" # different from the tutorial so I can have 2 versions 
sc.settings.set_figure_params(dpi = 80, frameon = False, figsize=(3,3), facecolor='white') # dpi = dots per inch

# low dpi = yields small inline figures

adata = sc.datasets.paul15() 
adata 

sc.pp.recipe_zheng17(adata)
# Normalization and filtering

# sc.tl.diffmap(adata)
# sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')

sc.tl.pca(adata, svd_solver = 'arpack')
# arpack is the wrapper for scipy 
# svd = single value decomposition for matrices
# factorization of a real or complex matrix 

sc.pp.neighbors(adata, n_neighbors=4, n_pcs = 20)
sc.tl.draw_graph(adata)

# sc.pl.draw_graph(adata, color="paul15_clusters", legend_loc='on data')
# ^^ a bit messy

# sc.pl.highest_expr_genes(adata, n_top=10)

sc.tl.leiden(adata, resolution = 1.0)
sc.tl.paga(adata, groups="leiden")
# sc.pl.paga(adata, color=["leiden", "Hba-a2", "Elane", "Irf8"])
# why did they decide these specifically? I tried any otehr combo and it doesn't work

# sc.pl.paga(adata, color=['leiden', 'Itga2b', 'Prss34', 'Cma1'])

adata.obs['leiden'].cat.categories
# Index(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20'],  dtype='object')

adata.obs['leiden_anno'] = adata.obs['leiden']
# print (adata.obs['leiden'])

adata.obs['leiden_anno'].cat.categories = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10/Ery', '11', '12','13', '14', '15', '16/Stem', '17', '18', '19/Neu', '20/Mk'] # apparently the extra 4 shouldn't be there

sc.tl.paga(adata, groups ='leiden_anno')
sc.pl.paga(adata, threshold=0.03, show=False)

sc.tl.draw_graph(adata, init_pos = 'paga')
# Notes for component documentation


There are `.read_[]x_[]` formats for converting to a data format that is understood

- Var_names (subset of `adata.var_names`) — key used as label to group the values (`var_group_labels`) and mapping values for sequences of valid `adata.var_names`...
- writing an h5ad cache file to speedup reading next time

## Preprocessing (`pp`)

### `pp.neighbors`
neighborhood graph of observations
relies on UMAP -- provies method for estimating connectivities of data points 

- `n_neighbors` = size of local nieghborhood 
- `n_pcs` = # of PC(principal components)s 

### `pp.calculate_qc_metrics`
- calculates qc metrics: `qc_vars` is a Collection[str] and stores the key for the columns of `.var` that you control for
- `inplace` = to place caluclated metrics in .obs or .var (??)
- `percent_top` = which proportion of top genes to cover (if none, don’t calc) — for example, `percent_top=[50]` finds cumulative proportion to the 50th most expressed gene
- `log1p` sets to false by skip computing log1p transformed annotations


```python
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
```

^ subset of `.calculate_qc_metrics` → the # of genes with at least 1 count in a cell, calc of all cells (`n_genes_by_{expr}`)


### `pp.scale`
```python
sc.pp.scale(adata, max_value=10)
```

scale data to unit variance and zero mean  ([https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.scale.html#scanpy-pp-scale](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.scale.html#scanpy-pp-scale))

### `.raw.to_data()` (not sure what cat)
the `.raw` attribute of the anndata set to *normalized* and *logarithmized* raw gene expression — freezes the state of AnnData object 

you can get back anndata of the object in `.raw` by calling `.raw.to_adata()` (this converst anncollection to annadata object)

- If you don’t proceed below with correcting the data with `sc.pp.repress_out` and scaling it via `sc.pp.scale`
    - dont need to use .raw
- result of the *highly varialve genes* is stored as an annotation in the `.var.highly_variable` and autodetected by PCA and therefore `sc.pp.neighbors` and any other mainfold/graph tools
    - `sc.pp.neighbors`: compute neighborhood graph of observations
    - neighbor search efficiency of this relies on UMAP → method for estimating *connectivities* of data points (`method==’umap’`) if `method=='gauss'` then computed according to a different algorithm ([[Coifman05]](https://scanpy.readthedocs.io/en/latest/references.html#coifman05) [[Haghverdi16]](https://scanpy.readthedocs.io/en/latest/references.html#haghverdi16))


### `.external.pp.bbknn`
Batch-balanced kNN -- alters the kNN procedure to identify each cell's top neighbors in each batch separately instead of the entire cell bool with no accounting for batch
used as an alternative to `neighbors()`. 

## Plotting (`pl`)
makes violin plots — wraps `seaborn.violinplot()` for AnnData 

- the kets are for accessing the variables of `.var_names` or fields of `.obs`
- jitter (add to stripllot only when its true)
- multi_panel: display the guys in multiple panels (different graphs)

## Tools (`tl`)

###  `tl.paga`
maps out the course-grained connectivity structures of complex manifolds (from documentation file:///Users/clairebookworm/Documents/scanpy/scanpy.readthedocs.io/en/latest/generated/scanpy.tl.paga.html)

PAGA = **partition-based graph abstraction** -> through quantifying the connectivity of partitions in the single-cell graph, it generals a much simpler abstracted graph (*paga graph*) of partitions -> edge weights rep confidence in the presence of connetions 

`paga()` -> thresholding this confidence here makes better *manifold* data (what's manifold data?)
confience is interpreted as a ratio of the actual vs. expected value of connections (no p-value)

- groups = key for categorical `adata.obs` and you can pass in predefined groups (default is first key of leiden or louvain
- use_rna_velocity = RNA velocity to orient edges (rna velocity?)
= model = version of paga connectivity model

### tl.leiden
(this is from the documentation, whatever the examples here had)
cluster cells into subgroups
requires `neighbors()` or `bbknn()` having run first 

- adjacency = sparse adjacency matrix of graph
- directed = is the graph directed
- use_weights = edges in graph are weighted
- neighbors_key = neighbor connectively as adjacency (any edge?) -> no specified looks at `.obsp['connectivities']`
- obsp = use obsp[obsp] as adjacency


### `tl.pca` 
is principal component analysis. it computes PCA coordinates, loadings, and variance decomposition. 
- `zero_center`: if true, compute standard PCA from covar matrix, if false, omit zero-centering variables (which allows to handle sparse input efficiently.) none is auto based on sparseness of data
- `svd_solver`: ARPACK (wrapper in SciPy), randomized (rando algo from halko 2009), auto (chooses auto depending on size of the problem), logpcg (alternative scipy solver)

(svd = single value denominator))

### `tl.ingest`
map labels and embeddings from ref data to new data 
integrates the embed and annotations in `adata` using the `adata_ref` dataset through projecting on a PCA that has been fitted on the reference data
Uses a knn calssified for mapping labels and UMAP for mapping embeddings 

`obs` = labels keys in `adata_ref.obs` which need to be mapped to `adata.obs`
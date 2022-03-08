# March 2

h5ad file

- hdf5 file with some additional structure specifying how to store AnnData objects

```python
# imports 
import numpy as np
import pandas as pd
import scanpy as sc
```

# **Notes**

```python
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
```

```python
adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading
```

There are `.read_[]x_[]` formats for converting to a data format that is understood

- Var_names (subset of `adata.var_names`) — key used as label to group the values (`var_group_labels`) and mapping values for sequences of valid `adata.var_names`...
- writing an h5ad cache file to speedup reading next time

```python
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`

adata # shows u the info about the adata (which is a sc.AnnData object from the mtx file)
```

AnnData object with n_obs × n_vars = 2700 × 32738
var: 'gene_ids'

![Screen Shot 2022-03-02 at 1.00.22 PM.png](March%202%20c3163/Screen_Shot_2022-03-02_at_1.00.22_PM.png)

there’s a lot of mess data here 

```python
#  basic filtering 
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
```

> filtered out 19024 genes that are detected in less than 3 cells
> 

[https://master.bioconductor.org/packages/release/workflows/html/simpleSingleCell.html#examining-gene-level-metrics](https://master.bioconductor.org/packages/release/workflows/html/simpleSingleCell.html#examining-gene-level-metrics)

```python
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
```

[https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.calculate_qc_metrics.html#scanpy-pp-calculate-qc-metrics](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.calculate_qc_metrics.html#scanpy-pp-calculate-qc-metrics) 

- calculates qc metrics: `qc_vars` is a Collection[str] and stores the key for the columns of `.var` that you control for
- `inplace` = to place caluclated metrics in .obs or .var (??)
- `percent_top` = which proportion of top genes to cover (if none, don’t calc) — for example, `percent_top=[50]` finds cumulative proportion to the 50th most expressed gene
- `log1p` sets to false by skip computing log1p transformed annotations

```python
import scanpy as sc
import seaborn as sns

pbmc = sc.datasets.pbmc3k()
pbmc.var["mito"] = pbmc.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(pbmc, qc_vars=["mito"], inplace=True)
sns.jointplot(
    data=pbmc.obs,
    x="log1p_total_counts",
    y="log1p_n_genes_by_counts",
    kind="hex",
)
```

Calculates a number of qc metrics for an AnnData object, see section Returns for specifics. Largely based on calculateQCMetrics from scater [McCarthy17]. Currently is most efficient on a sparse CSR or dense matrix.

A violin plot of some of the computed quality measures:

- the number of genes expressed in the count matrix
- the total counts per cell
- the percentage of counts in mitochondrial genes

```python
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
```

makes violin plots — wraps `seaborn.violinplot()` for AnnData 

- the kets are for accessing the variables of `.var_names` or fields of `.obs`
- jitter (add to stripllot only when its true)
- multi_panel: display the guys in multiple panels (different graphs)

```python
import scanpy as sc
adata = sc.datasets.pbmc68k_reduced()
sc.pl.violin(adata, keys='S_score')
# another example
sc.pl.violin(adata, keys='S_score', groupby='bulk_labels', rotation=90)
# another
groupby_order = ['CD34+', 'CD19+ B']
sc.pl.violin(adata, keys='S_score', groupby='bulk_labels', rotation=90,
    order=groupby_order)
# hmm
sc.pl.violin(adata, keys='S_score', stripplot=False)
```

stripplot false means no dots around

![Untitled](March%202%20c3163/Untitled.png)

![Untitled](March%202%20c3163/Untitled%201.png)

![Untitled](March%202%20c3163/Untitled%202.png)

```python
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
```

obviously scatter plot , putting in `adata` and the x y axes

```python
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
```

^ subset of `.calcualte_qc_metrics` → the # of genes with at least 1 count in a cell, calc of all cells (`n_genes_by_{expr}`)

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata) # log's the data 

# OUTPUT 
# normalizing counts per cell
# /Users/clairebookworm/opt/miniconda3/lib/python3.9/site-packages/scanpy/preprocessing/_normalization.py:155: UserWarning: Revieved a view of an AnnData. Making a copy.
#   view_to_actual(adata)
#     finished (0:00:00)
```

**highly-variable** genes 

- plot dispersions or normalized variance vs means for genes

```python
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# extracting highly variable genes
#     finished (0:00:02)
# --> added
#     'highly_variable', boolean vector (adata.var)
#     'means', float vector (adata.var)
#     'dispersions', float vector (adata.var)
#     'dispersions_norm', float vector (adata.var)

sc.pl.highly_variable_genes(adata)
adata.raw = data 
```

^ not entirely sure what this means, is this from the paper: [https://www.nature.com/articles/ncomms14049](https://www.nature.com/articles/ncomms14049) 

meanvarplot and variablefeatureplot

![Untitled](March%202%20c3163/Untitled%203.png)

the `.raw` attribute of the anndata set to *normalized* and *logarithmized* raw gene expression — freezes the state of AnnData object 

you can get back anndata of the object in `.raw` by calling `.raw.to_adata()` (this converst anncollection to annadata object)

**pp = preprocessing!** 

- If you don’t proceed below with correcting the data with `sc.pp.repress_out` and scaling it via `sc.pp.scale`
    - dont need to use .raw
- result of the *highly varialve genes* is stored as an annotation in the `.var.highly_variable` and autodetected by PCA and therefore `sc.pp.neighbors` and any other mainfold/graph tools
    - `sc.pp.neighbors`: compute neighborhood graph of observations
    - neighbor search efficiency of this relies on UMAP → method for estimating *connectivities* of data points (`method==’umap’`) if `method=='gauss'` then computed according to a different algorithm ([[Coifman05]](https://scanpy.readthedocs.io/en/latest/references.html#coifman05) [[Haghverdi16]](https://scanpy.readthedocs.io/en/latest/references.html#haghverdi16))

```python
adata = adata[:, adata.var.highly_variable]
 # regress out effects of TOTAL 
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# regressing out ['total_counts', 'pct_counts_mt']
#     sparse input is densified and may lead to high memory use
# /Users/clairebookworm/opt/miniconda3/lib/python3.9/site-packages/statsmodels/tsa/base/tsa_model.py:7: FutureWarning: pandas.Int64Index is deprecated and will be removed from pandas in a future version. Use pandas.Index with the appropriate dtype instead.
#   from pandas import (to_datetime, Int64Index, DatetimeIndex, Period,
# /Users/clairebookworm/opt/miniconda3/lib/python3.9/site-packages/statsmodels/tsa/base/tsa_model.py:7: FutureWarning: pandas.Float64Index is deprecated and will be removed from pandas in a future version. Use pandas.Index with the appropriate dtype instead.
#   from pandas import (to_datetime, Int64Index, DatetimeIndex, Period,
#     finished (0:00:15)
```

then we scale the gene to unit **variance** - clip values exceeding stdev 10

```python
sc.pp.scale(adata, max_value=10)
```

scale data to unit variance and zero mean  ([https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.scale.html#scanpy-pp-scale](https://scanpy.readthedocs.io/en/latest/generated/scanpy.pp.scale.html#scanpy-pp-scale)) 

### Principal component analysis

# **Trial and Error**

Sc.write_10x_h5

```jsx
➜  data python3 -c "import scanpy as sc"
Traceback (most recent call last):
  File "<string>", line 1, in <module>
  File "/usr/local/lib/python3.9/site-packages/scanpy/__init__.py", line 19, in <module>
    from anndata import AnnData, concat
ImportError: cannot import name 'concat' from 'anndata' (/usr/local/lib/python3.9/site-packages/anndata/__init__.py)

➜  data pip install -e .
pip install "anndata<=0.7.3"
python3 -c "import scanpy as sc"

➜  data pip3 install scanpy   --target $PYTHONPATH --upgrade
```

Attempts at figuring out how to read a h5ad file:

```python
import numpy as np
import pandas as pd
import scanpy as sc
 
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
 
sc.read_10x_h5(
   'sciPlex1_HEK293T.h5ad',
)
```

As well as tutorial clustering code:

```python
import matplotlib
matplotlib.use('Agg')  # plotting backend compatible with screen
import sys
import scanpy as sc
 
sc.settings.verbosity = 2  # show logging output
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=300)  # set sufficiently high resolution for saving
 
filename = sys.argv[1]  # read filename from command line
 
def basic_analysis(filename):
   adata = sc.read_10x_h5(filename)
   sc.pp.recipe_zheng17(adata)
   sc.pp.neighbors(adata)
   sc.tl.louvain(adata)
   sc.tl.paga(adata)
   sc.tl.umap(adata)
   sc.tl.rank_genes_groups(adata, 'louvain')
   adata.write('./write/result.h5ad')
   # plotting
   sc.pl.paga(adata)
   sc.pl.umap(adata, color='louvain')
   sc.pl.rank_genes_groups(adata, save='.pdf')
 
 
if __name__ == "__main__":
   basic_analysis(filename)
```

*Still trying to figure out everything*… (2/22)

# **Data File Formats**

[https://broadinstitute.github.io/wot/file_formats/](https://broadinstitute.github.io/wot/file_formats/)

## anndata

![Screen Shot 2022-03-02 at 1.08.38 PM.png](March%202%20c3163/Screen_Shot_2022-03-02_at_1.08.38_PM.png)

`obsp` is **pairwise annotations of observations**, a mutable mapping with array-like values

- stores for each key 2+ dimensional `ndarray` where first 2 dimesions are of length `n_obs` — sliced with `data` and `obs` but behaves like a mapping

`anndata.AnnData.var`

- **1d annotation of vars and features**
- returns `pd.DataFrame`
    - dataframe is in `pandas`, a data structure with labeled axes and arithmetic operations on both row and column

`anndata.AnnData.obsm` — **multi-dimensional annotation** of observations 

- stores for each key 2+ `ndarray` and is v similar to obsp

## **MTX (MetaStream ASCII text format)**

%%MatrixMarket matrix coordinate real general%32738 2700 228688432709 1 432707 1 132706 1 1032704 1 132703 1 532702 1 6

- uses XML-style tags to store information about the scene, such as
    - global scene options
    - object (instance) hierarchies
    - camera locations
    - animation data

## ***h5ad / hdf5***

A HDF5 file that provides a scalable way of keeping track of data together with learned annotations.

[https://anndata.readthedocs.io/en/latest/](https://anndata.readthedocs.io/en/latest/)

**Sample**: [https://www.covid19cellatlas.org/](https://www.covid19cellatlas.org/) (h5ad file format for covid)

âHDFˇˇˇˇˇˇˇˇL∫ˇˇˇˇˇˇˇˇ`à(r)à(r)TREEˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇ‡HEAPX »Xobsvar8(-TREEˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇˇËèHEAPX àdataindicesindptr8SNOD HhPìxìòïg#yt|yî~yàHencoding-typePencoding-versionGCOLcsc_matrix0.1.0	dataframe0.1.0pvalqvaltop_to_second_best_ratio	top_oligo		cell_typen.umi	hash_umisSize_FactorsampleCell_indexE03_C05_C05E03_C05_C09E03_C05_D05

## ***Loom***

Hdf5 file format for efficient storage.
# Scanpy Session 1 Notes

bulk = smoothie (same as scs analogy)

**examples of scs sequencing methods**: 
- CEL-seq2/C1
- Drop-seq
- MARS-seq
- SCRB-seq
- Smart-seq/C1
- Smart-seq2

there are some full-lenngth (whole transcript: smart-seq) and 3'biased protocol (captures 3' end, UMI): the rest 

UMI = unique molecular identifier 

Libraries: Feature, GEX
- sequencing: BCL to mkfastq
- cell ranger that comes with the 10x genomics system that helps you map your data
- FASTQ file 
	- first line that starts with the @
	- followed by sequence id and optional descip of file
	- line 2 is raw sequence letter read out in the ngscore facility
	- line 3 starts with plus and has some description of the data

## Sensitivity of single-cell RNAseq 
molecular detection limit -> how many copies do I need to detect the gene at least once 

smartseq is sensitive vs. dropseek, chromium, 10x genomics that are high thruput but more shallow 

## Workflow 
now we have our **count matrix** = dimension cells x genes (30k genes, 3k genes in the example)
- every entry is the # of RNA molecules detected 
- huge data set but you don't want to go thru one by one 

PREPROCESSING (80% of time)
1. **QC** distinguish if the cell are actualyl cells or debris
2. normalizing cells for cell-specific and cell-type specific bias and gene-length (if full lenght protocol)
3. feature select (like highly variable genes) - if the variability of the gene is high then the intrinsice noise of the experience is fixed
4. data correction (correction for systemic bias across samples)

DIMENSIONALITY REDUCTION and VIS
- 30k dimensions to 2 -> 2dimensional

*downstream analysis*
1. clustering (group the data by dissimilar properties)
2. annotate (using marker genes - prior knolwedge to see if a specific gene is expression to assign name to a cluster)
3. trajectory inference (how the stem cells change to become the respective cell types -- where do they undergo decisions and what genes are involved, pseudo-time)
4. gene dynamics 
5. differential expression (what is changing in cell upon disease vs healthy - which genes change across groups)
6. compositional analysis (how do compositions of my cells change in different conditions -- capture if there are cell types with different appearance)

**(Luecken, Theis paper on Current best practices)**

## Scanpy
uses Anndata 
- `.var` = anno of variables / features (pandas dataframe)
- `.x` = data matrix of shape - #obserations x variables (numpy array, scipy sparse matrix)
- `.obs` = anno of obseervations (pandas dataframe)
- `.uns` = unstructured anno (dict) 
	- anything that doesn't fit like colors of a plot


```
AnnData object with n obs x n vars 1004 x 12220 
	obs: 'n counts', 'log_ counts', 'n_genes', 'mt_frac', 'size_factors', 'S_score', 'G2M score', 'phase', 'louvain_r 1.5', 'louvain_r0.5', 'annotated' 
	var: 'gene ids', 'n_ cells', 'highly_variable', 'means' , 'dispersions', 'dispersions norm' 
	uns: 'pca', 'neighbors', 'diffmap_evals', 'draw _graph', 'phase_colors', 'louvain', 'louvain_r1.5_colors', 'louvai n_r0.5_colors', 'rank_genes_r0.5', 'rank _genes_r1.5', 'annotated _colors', 'rank_genes_groups', paga', annotated siz es 
	obsm: 'X pa', 'X tone', 'X umap', 'X diffmap', 'X draw graph fa' 
	varm: 'PCs 
	layers: 'counts'
```

^ example of AnnData

`obsm` contains other dimensions of objects; `varm` components; `layers` same dimensions of adata object, count of every molecule (x matrix might have different)

downsampling to: `sc.pp.subsample(adata,n_obs=5000)`

```
mkdir /broad/hptmp/cwang
curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
export PATH=$PATH:/broad/hptmp/cwang
export PATH=$PATH:/home/unix/clwang
export PATH=$PATH:/home/unix/clwang/.local/bin
python get-pip.py
pip install numpy 
```
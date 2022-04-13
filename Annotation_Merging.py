# import statements
import collections
import numpy as np
import pegasus as pg
import pandas as pd
import pegasusio as io
# cd -q ..
# cd -q -

# 2d arrayclear
    # Each row correspond to a dataset
    # Within that row, entry 1 is the path to the dataset, entry 2 is the path to the annotations, entry 3 is what the file should be named
paths = [
#     ["../../data/human/HTA/Epityphlon/GSM3980132_Adult-Epityphlon1_dge.txt", "../../output_pg/human_tissue_atlas/Duodenum/1.4-mad-2/", "human-HTA-Epityphlon"],
#         ["../../data/human/HTA/Ureter/GSM4008665_Adult-Ureter1_dge.txt", "../../output_pg/human_tissue_atlas/Ureter/1.4-mad-2/", "human-HTA-Ureter"],
        ["../../data/human/HTA/Uterus/GSM4008666_Adult-Uterus1_dge.txt", "../../output_pg/human_tissue_atlas/Uterus/1.4-mad-2/", "human-HTA-Uterus"]]

for path in paths:
    print("STARTING: " + path[2])
    data = pg.read_input(path[0]) # read data

    #read annotation files
    cells = pd.read_csv(path[1] + "!cells.csv")
    clusters = pd.read_csv(path[1] + "!clusters.csv")
    if 'Unnamed: 0' in cells.columns:
        del cells['Unnamed: 0']
    if 'Unnamed: 0' in clusters.columns:
        del clusters['Unnamed: 0']

    # formats barcodes in annotation to be in the same format as how they are in the actual data
    barcodes_proper = []
    for i, j in cells.iterrows():
        orig_barcode = cells["barcodekey"][i]
        channel = cells["Channel"][i] + "-"
        barcodes_proper.append(orig_barcode.replace(channel, ""))
    cells["barcodes_proper"] = barcodes_proper

    # currently, the barcodes in data.obs are the indexes. It is easier to access them if they are a separate column. This accomplishes that.
    data_barcodes = []
    for barcode in data.obs.index:
        data_barcodes.append(barcode)
    data_barcodes
    data.obs['barcodes'] = data_barcodes

    # Getting the cell types from the clusters df and putting it into the cells df
    celltype1 = []
    celltype2 = []
    for i, j in cells.iterrows():
        cluster = cells["louvain_labels"][i]
        cell1 = clusters["cell_type"][cluster-1]
        cell2 = clusters["cell_type2"][cluster-1]
        celltype1.append(cell1)
        celltype2.append(cell2)
    cells["cell_type"] = celltype1
    cells["cell_type2"] = celltype2

    # subsetting out cells in the data that are not present in the annotations
    subset = io.MultimodalData(data[data.obs['barcodes'].isin(barcodes_proper)].copy())

    # sorting annotations based on ordering of data.obs
    data_barcodes = data.obs['barcodes'].tolist()
    cells.barcodes_proper = cells.barcodes_proper.astype("category")
    cells.barcodes_proper.cat.set_categories(data_barcodes, inplace=True)
    cells = cells.sort_values(["barcodes_proper"])

    # getting relevant annotation columns into data.obs
    subset.obs['cluster'] = cells['louvain_labels'].values.astype(str)
    subset.obs['cell_type'] = cells['cell_type'].values
    subset.obs['cell_type2'] = cells['cell_type2'].values
    subset.obs['pca1'] = cells['pca1'].values.astype(str)
    subset.obs['pca2'] = cells['pca2'].values.astype(str)
    subset.obs['umap1'] = cells['umap1'].values.astype(str)
    subset.obs['umap2'] = cells['umap2'].values.astype(str)

    pg.write_output(subset, path[2] + ".h5ad")
    print("COMPLETED: " + path[2])
    print()
    print()
    print()

# scp clwang@login.broadinstitue.org:/broad/hptmp/cwang/ipynb /

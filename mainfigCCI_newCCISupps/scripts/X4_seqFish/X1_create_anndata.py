"""
Creating anndata object from the seqFISH data.

    INPUT: * /Volumes/GML001-Q1851/Brad/seqFISH/cortex_svz_*.csv
    OUTPUT: * /Volumes/GML001-Q1851/Brad/seqFISH/svz.h5ad
"""

################################################################################
                     # Environment Setup #
################################################################################
# TODO change this to your directory
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/stlearn_manuscript/mainfigCCI_newCCISupps/'

import os
os.chdir(work_dir)

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from PIL import Image

data_dir = '/Volumes/GML001-Q1851/Brad/seqFISH/'
out_dir = data_dir

################################################################################
                    # Loading in the data #
################################################################################
counts = pd.read_csv(data_dir+'cortex_svz_counts.csv', sep=',')
spatial = pd.read_csv(data_dir+'cortex_svz_cellcentroids.csv', sep=',')
spatial.index = spatial.index.values.astype(str)
cell_annots = pd.read_csv(data_dir+'cortex_svz_cell_type_annotations.csv',
                                                           sep=',', index_col=0)
cell_annots.index = cell_annots.index.values.astype(str)

################################################################################
                    # Creating the anndata #
################################################################################
data = ad.AnnData(counts)

data.obs["imagecol"] = -spatial["Y"].values #(max(spatial["Y"])-spatial["Y"]).values
data.obs["imagerow"] = spatial["X"].values

# Create image
quality = 'hires'
max_size = np.max([data.obs["imagecol"].max(),data.obs["imagerow"].max()])
max_size = int(max_size + 0.1*max_size)

image = Image.new('RGB', (max_size, max_size), (255, 255, 255))
imgarr = np.array(image)

library_id = "seqFISH+"

# Adding image to anndata #
data.uns["spatial"] = {}
data.uns["spatial"][library_id] = {}
data.uns["spatial"][library_id]["images"] = {}
data.uns["spatial"][library_id]["images"][quality] = imgarr
data.uns["spatial"][library_id]["use_quality"] = quality
data.uns["spatial"][library_id]["scalefactors"] = {}
data.uns["spatial"][library_id]["scalefactors"]\
                                           ["tissue_" + quality + "_scalef"] = 1
data.uns["spatial"][library_id]['scalefactors']['spot_diameter_fullres'] = 60
data.obsm["spatial"] = spatial[["Y","X"]].values
data.obsm["spatial"][:,0] = -data.obsm["spatial"][:,0]

# Adding in cell annotation information #
data.obs['cluster'] = cell_annots.loc[:,'louvain'].astype('category')
data.obs['region'] = spatial.loc[:,'Region']
data.obs['region'] = data.obs['region'].astype('category')

# Subsetting to SVZ #
svz = data[data.obs['region'].values=='SVZ',:].copy()
sc.pl.spatial(svz, color='cluster')

### Adding the cell type information in based on paper ###
ct_map = {'1': 'Choroid Plexus', '14': 'Oligodendrocytes',
          '20': 'Microglia', '7': 'Interneuron',
          '10': 'Astrocytes', '2': 'Endothelial',
          '22': 'Ependymal', '8': 'neural-stem',
          '16': 'neural-stem', '15': 'TAPs'}
cell_types = [f'unknown_{clust}' if clust not in ct_map else ct_map[clust]
              for clust in svz.obs['cluster'].values.astype(str)]
svz.obs['cell_type'] = cell_types
svz.obs['cell_type'] = svz.obs['cell_type'].astype('category')

sc.pl.spatial(svz, color='cell_type')

################################################################################
                    # Gene filtering & normalisation #
################################################################################
sc.pp.filter_genes(svz, min_cells=3)
sc.pp.normalize_total(svz)

################################################################################
                    # Saving the data #
################################################################################
#svz.write_h5ad(out_dir+'svz.h5ad', compression='gzip')
svz.write_h5ad(out_dir+'svz2.h5ad', compression='gzip')




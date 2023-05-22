"""
 Loads in the slide-seq data & subsequently formats the spatial
                information and normalises data using scanpy workflow prior to
                downstream lr-cci analysis.

                 INPUT: * /Volumes/GML001-Q1851/Brad/slideSeq/hipp.h5ad
                 OUTPUT: * /Volumes/GML001-Q1851/Brad/slideSeq/hipp.h5ad
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

data_dir = '/Volumes/GML001-Q1851/Brad/slideSeq/'
out_dir = data_dir

################################################################################
                        # Loading the data #
################################################################################
data = sc.read_h5ad(data_dir+'hipp.h5ad')

################################################################################
                     # Performing normalisation #
################################################################################
print(data.shape)
sc.pp.filter_genes(data, min_cells=3)
sc.pp.filter_cells(data, min_counts=50)
print(data.shape)
sc.pp.normalize_total(data)

################################################################################
                 # Reformattting spatial information #
################################################################################
# Create image
quality = 'hires'
max_size = np.max([data.obs["imagecol"].max(),data.obs["imagerow"].max()])
max_size = int(max_size + 0.1*max_size)

image = Image.new('RGB', (max_size, max_size), (255, 255, 255))
imgarr = np.array(image)

library_id = "slideSeq"

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
data.obsm["spatial"] = data.obs.loc[:,["imagecol","imagerow"]].values

# Plotting to check #
sc.pl.spatial(data, color='TTR', size=.8)
sc.pl.spatial(data, color='cell_type', size=.8)

# Looks great !!! #
#data.write_h5ad(out_dir+'hipp.h5ad', compression='gzip')
data.write_h5ad(out_dir+'hipp_rep.h5ad', compression='gzip')








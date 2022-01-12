"""
Performs the preprocessing on the skin data.

                INPUT: * /Volumes/GML001-Q1851/Jon/C1_skin_raw/
                OUTPUT: * /Volumes/GML001-Q1851/Brad/skin_analysis/
"""

################################################################################
                        # Environment setup  #
################################################################################
# TODO: NOTE must be run in folder with README entitled StLearn Reproduce
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/'
import os
os.chdir(work_dir)

import numpy as np
import pandas as pd
import scanpy as sc
import stlearn as st
import matplotlib.pyplot as plt

import scripts.utils.helpers as g_hs

gml_rdm_path = '/Volumes/GML001-Q1851/'
out_dir = '/Volumes/GML001-Q1851/Brad/skin_analysis/'

################################################################################
                        # Load the data  #
################################################################################
#### Loading the skin data so can perform the significance testing #####
data = g_hs.load_skin()

##### Subsetting the skin data to just the spots from one sample #######
data.obsm['xs'] = data.obsm['spatial'][:,0]
data.obsm['ys'] = data.obsm['spatial'][:,1]

st.pl.het_plot(data, use_het='xs')
plt.show()

st.pl.het_plot(data, use_het='ys')
plt.show()

slice = np.array([0]*data.shape[0])
bool_ = data.obsm['xs'] > 850
bool_2 = np.logical_and(data.obsm['xs'] > 650,
                        data.obsm['ys'] > 1800)
bool_3 = np.logical_and(data.obsm['xs'] > 800,
                        data.obsm['ys'] > 1750)
bool_ = np.logical_or(bool_, bool_2)
bool_ = np.logical_or(bool_, bool_3)
slice[bool_] = 1

slice = np.array([f'slice {i}' for i in slice])

data.obs['slice'] = pd.Series(slice, index=data.obs_names).astype('category')

# Subsetting to just the first slice ! #
slice = data.obs['slice'].values
data = data[slice=='slice 0',:].copy()

# After the spot subsetting, remove genes which aren't expressed #
print(data.shape)
st.pp.filter_genes(data, min_cells=3)
print(data.shape)

sc.pl.spatial(data, color='MITF', cmap='Spectral_r')

#################### Saving the data #############################
data.write_h5ad(out_dir+'skin_ppd.h5ad', compression='gzip')

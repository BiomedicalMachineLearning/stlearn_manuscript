"""
Adds in the stLearn clustering information to the anndata generate from SME
clustering.

         INPUT:  * /Volumes/GML001-Q1851/Brad/
                                   breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad
                 * data/breast/bc_clusters.csv
         OUTPUT: * /Volumes/GML001-Q1851/Brad/breast_ClusterCCIResults.h5ad
                 * plots/X2_method_comp_clusters/
"""

################################################################################
                    # Environment setup #
################################################################################
#TODO: NOTE must be run in folder with README entitled StLearn Reproduce
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/'
stlearn_path = '/Users/uqbbalde/Desktop/Uni_Studies/myPython/' \
               'stlearn_latest/stLearn/'
gml_rdm_path = '/Volumes/GML001-Q1851/'
data_dir = gml_rdm_path+'Brad/'

import os, sys
os.chdir(work_dir)
sys.path.append(stlearn_path)

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import stlearn as st
st.settings.set_figure_params(dpi=180)
from pathlib import Path

import stlearn.plotting.cci_plot as cci_plot

import beautifulcells.visualisation.helpers as vhs
import beautifulcells.visualisation.quick_plots as qpl
qpl.setUp()

mode = '' # 'abs_' is discrete spot labels, '' is deconvolution.
data_set = 'breast'
data_dir2 = 'data/breast/'
out_plots = 'plots/X2_method_comp_clusters/'

################################################################################
 # Loading data, performing clustering, & calculate cell type interactions ! #
################################################################################
# Loading the data #
data = sc.read_h5ad(data_dir + f'breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad')

###### Adding in the cluster labels #####
bc_meta = pd.read_csv(data_dir2+'bc_clusters.csv', index_col=0)
overlap_bcs = [bc for bc in bc_meta.index if bc in data.obs_names]
bc_meta = bc_meta.loc[overlap_bcs,:]
data = data[overlap_bcs,:]
data.obs['cell_type'] = bc_meta.loc[:,'louvain'].values.astype(str)

##### Plotting & setting colors ######
sc.pl.spatial(data, color='cell_type', show=False)
vhs.dealWithPlot(True, True, True, out_plots, 'breast_clusters_spatial.pdf',
                 300)
# Setting colors to be equivalent to cellchat colors #
colors = {'0': '1F77B4', '1': 'F87F13', '2': '359C62',
          '3': 'D32929', '4': '69308E',
          '5': '8C564C', '6': 'F33CA9',
          '7': 'F5E801', '8': '08F7F0', '9': 'AFC7E6', '10': 'FEBB7A'}
label_set = list(data.obs['cell_type'].cat.categories)
label_colors = ['#'+colors[label] for label in label_set]
data.uns['cell_type_colors'] = label_colors

sc.pl.spatial(data, color='cell_type', show=False)
vhs.dealWithPlot(True, True, True, out_plots,
                 'breast_clusters_spatial.pdf', 300)

# Counting the LR interactions per cell type #
st.tl.cci.run_cci(data, 'cell_type', min_spots=3,
                  spot_mixtures=False, n_perms=1000) # NOTE not significant in mixture mode.
print(data.uns['lr_cci_raw_cell_type'].values.sum())

# Writing this out #
data.write_h5ad(data_dir+'breast_ClusterCCIResults.h5ad', compression='gzip')




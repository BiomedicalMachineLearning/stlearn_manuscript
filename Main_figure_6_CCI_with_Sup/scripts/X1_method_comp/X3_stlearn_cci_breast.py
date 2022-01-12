"""
(21/09/21) -> Creating the equivalent plots as CellChat for stlearn, in order to
                do a side-by-side comparison, which if considered in conjunction
                with the colocalisation results & the spatial plots of location
                of the interactions, should show how stlearn can reduce the
                false-positive-rate by taking into account spatial information.

                INPUT:  * /Volumes/GML001-Q1851/Brad/
                               breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad
                        * data/label_transfer_bc.csv
                OUTPUT: * plots/method_comp/stlearn_breast_*
"""

################################################################################
                # Environment setup & load test data #
################################################################################
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

import stlearn.plotting.cci_plot as cci_plot

import scripts.utils.visualisation.helpers as vhs
import scripts.utils.visualisation.quick_plots as qpl
qpl.setUp()

mode = '' # 'abs_' is discrete spot labels, '' is deconvolution.
data_set = 'breast'
out_plots1 = 'plots/method_comp/'
out_plots = 'plots/method_comp/stlearn_breast_'

################################################################################
 # Loading data, adding label transfer, & calculate cell type interactions ! #
################################################################################
# Loading the data #
data = sc.read_h5ad(data_dir+
                    f'breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad')
data.obsm['spot_neighbours'] = data.uns['spot_neighbours']

# Adding the label transfer information #
spot_mixtures = pd.read_csv('data/label_transfer_bc.csv', index_col=0, sep='\t')
labels = spot_mixtures.loc[:,'predicted.id'].values.astype(str)
spot_mixtures = spot_mixtures.drop(['predicted.id','prediction.score.max'],
                                   axis=1)
spot_mixtures.columns = [col.replace('prediction.score.', '')
                         for col in spot_mixtures.columns]
data.obs['cell_type'] = labels

data.uns['cell_type'] = spot_mixtures

# Plotting #
sc.pl.spatial(data, color='cell_type')
# Setting colors to be equivalent to cellchat colors #
colors = {'Bcell': '3FA239', 'Tcell': '469ED7', 'basal_like_1': 'D90C19',
          'basal_like_2': '246AA9', 'endothelial': '874293',
          'luminal_ar': 'E37C1A', 'macrophage': 'E46BA5',
          'mesenchymal': 'AD8ABE', 'stroma': '94431E'}
label_set = list(data.obs['cell_type'].cat.categories)
label_colors = ['#'+colors[label] for label in label_set]
data.uns['cell_type_colors'] = label_colors

sc.pl.spatial(data, color='cell_type', show=False)
vhs.dealWithPlot(True, True, True, out_plots1,
                 'breast_celltypes_spatial.pdf', 300)

# Counting the LR interactions per cell type #
data_sub = data[0:50,:]
st.tl.cci.run_cci(data, 'cell_type', min_spots=3,
                  spot_mixtures=True, n_perms=1000) # NOTE not significant in mixture mode.

print(data.uns['lr_cci_raw_cell_type'].values.sum())

# Writing this out #
data.write_h5ad(data_dir+'breast_LR&CCIResults.h5ad',
                compression='gzip')

################################################################################
            # Calculating the cell type interactions ! #
################################################################################
out='.png'
for sig in [True, False]:
    lrs = ['TGM2_ADGRG1', 'GPC3_IGF1R', 'SPP1_ITGB1']

    # Network plot #
    pos_1 = cci_plot.ccinet_plot(data, 'cell_type', min_counts=0, font_size=20,
                                 return_pos=True, figsize=(11,10),
                                 sig_interactions=sig)
    vhs.dealWithPlot(True, True, True, out_plots,
                                              f'{sig}-overall_network{out}', 300)
    for best_lr in lrs:
        cci_plot.ccinet_plot(data, 'cell_type', best_lr, min_counts=0,
                             figsize=(11,10), pos=pos_1, font_size=20,
                             sig_interactions=sig, pad=.25
                          )
        vhs.dealWithPlot(True, True, True, out_plots,
                         f'{sig}-{best_lr}_network{out}', 300)

    # Chord diagram #
    cci_plot.lr_chord_plot(data, 'cell_type', show=False, sig_interactions=sig)
    vhs.dealWithPlot(True, True, True, out_plots,
                                                f'{sig}-overall_chord{out}', 300)

    lrs = ['TGM2_ADGRG1', 'GPC3_IGF1R', 'SPP1_ITGB1']
    for lr in lrs:
        cci_plot.lr_chord_plot(data, 'cell_type', lr, show=False,
                               sig_interactions=sig)
        vhs.dealWithPlot(True, True, True, out_plots,
                                                   f'{sig}-{lr}_chord{out}', 300)

# TODO put the above into a powerpoint #

# This was used in chord_debug.py #
#lr_int_dfs[lrs[0]].to_csv(f'data/test/{lrs[0]}_int_df.txt', sep='\t')

""" DONE! Looks good now. 
Sort the heatmaps of the number of interactions. 
"""




















""" Junk code
"""
##### Attempt with the chord package ######
# Ref: https://pypi.org/project/chord/
from chord import Chord
import matplotlib.image as img

from PIL import Image

# reading image
img = Image.open('plots/ccis.png')

img.thumbnail((30, 30), Image.ANTIALIAS)

# bicubic used for interpolation
imgplot = plt.imshow(img)
plt.show()

matrix = [
    [0, 5, 6, 4, 7, 4],
    [5, 0, 5, 4, 6, 5],
    [6, 5, 0, 4, 5, 5],
    [4, 4, 4, 0, 5, 5],
    [7, 6, 5, 5, 0, 4],
    [4, 5, 5, 5, 4, 0],
]

names = ["Action", "Adventure", "Comedy", "Drama", "Fantasy", "Thriller"]

png_str = Chord(matrix, names).render_png()

Chord(matrix, names).to_png('plots/temp.png')

img = Image.open('plots/temp.png')

plt.imshow(img)
plt.show()

image = Image.frombytes("RGB", (200, 200), png_str, "raw", "BGRX")

plt.imshow(image)
plt.show()

#### Attempt with the mne package #####
from mne.viz import plot_connectivity_circle

N = 20  # Number of nodes
node_names = [f"N{i}" for i in range(N)]  # List of labels [N]

# Random connectivity
ran = np.random.rand(N,N)
con = np.where(ran > 0.9, ran, np.nan)

fig, axes = plot_connectivity_circle(con, node_names,
    colormap='Blues', facecolor='white', textcolor='black', colorbar=False,
    linewidth=10)










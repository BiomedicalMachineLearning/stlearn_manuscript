"""
Runs stLearn CCI analysis on the slide-seq data after gridding.

INPUT: * /QRISdata/Q1851/Brad/slideSeq/hipp_LR-CCI_clustered.h5ad

OUTPUT: * /QRISdata/Q1851/Brad/slideSeq/hipp_gridded_45-45_LR-CCI.h5ad
        * plots/X5_slideSeq/
"""

# TODO change this to your directory
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/stlearn_manuscript/mainfigCCI_newCCISupps/'

import os
os.chdir(work_dir)

import numpy as np
import scanpy as sc
from timeit import time

import scripts.utils.visualisation.helpers as vhs

import stlearn as st

data_dir = '/QRISdata/Q1851/Brad/slideSeq/'
out_dir = data_dir
out_plots = 'plots/X5_slideSeq/'

###########################################################################
			# Loading the data #
###########################################################################
data = sc.read_h5ad(data_dir+'hipp_LR-CCI_clustered.h5ad')

###########################################################################
                        # Gridding the data #
###########################################################################
n_row, n_col=45, 45
st.pl.grid_plot(data, use_label='leiden', n_row=n_row, n_col=n_col,
                show=False, size=.1)
vhs.dealWithPlot(True, True, True, out_plots, 'gridded_clusters.pdf', 300)

grid_data = st.tl.cci.grid(data, n_row=n_row, n_col=n_col, use_label='leiden')

sc.pl.spatial(grid_data, color='n_cells', size=2)
vhs.dealWithPlot(True, True, True, out_plots, 'gridded_cell-counts.pdf', 300)

### Visualising gridded cell types ###
cluster_set = grid_data.obs['leiden'].cat.categories.values
cluster_set2 = data.obs['leiden'].cat.categories.values
colors2 = np.array( data.uns['leiden_colors'] )
colors = [colors2[cluster_set2==cluster][0] for cluster in cluster_set]
grid_data.uns['leiden_colors'] = colors
st.pl.cluster_plot(grid_data, use_label='leiden', size=50)
vhs.dealWithPlot(True, True, True, out_plots, 'gridded_leiden.pdf', 300)

grid_data.obsm['deconvolution'] = grid_data.uns['leiden']
st.pl.deconvolution_plot(grid_data, show_donut=False, spot_size=60,
                         colors=grid_data.uns['leiden_colors'], celltype_threshold=.2,
                        )
vhs.dealWithPlot(True, True, True, out_plots, 'gridded_leiden-mixtures.pdf', 300)
###########################################################################
                       # Running the LR analysis #
###########################################################################
# Using same LR tested on single cell data
lrs = data.uns['lr_summary'].index.values.astype(str)

# Running the analysis #
n_cpus=2
start = time.time()
st.tl.cci.run(grid_data, lrs,
                  min_spots = 5, #Filter out any LR pairs with no scores for less than min_spots
                  distance=0, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=10000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=n_cpus, # Number of CPUs for parallel. If None, detects & use all available.
                  )
end = time.time()
diff = end-start
print(f'{diff} secs')
print(f'{diff/60} mins')

file_ = open(out_plots+f'grid_col-{n_col}_row-{n_row}_time.txt', 'w')
file_.write(f'{diff/60} mins/{n_cpus} cpus')
file_.close()

# grid_data.write_h5ad(out_dir+f'hipp_gridded_{n_col}-{n_row}_LR-CCI.h5ad', compression='gzip')
grid_data.write_h5ad(out_dir+f'hipp_gridded_{n_col}-{n_row}_LR-CCI_rep.h5ad',
                     compression='gzip')


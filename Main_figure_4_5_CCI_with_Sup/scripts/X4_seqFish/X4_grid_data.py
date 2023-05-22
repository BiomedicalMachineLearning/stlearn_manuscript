"""
Runs gridding function on the seqFISH+ data, essentially pseudobulking by
 neighbourhoods to improve speed with running CCI.
 Then runs stlearn LR-CCI analysis.

    INPUT: * /Volumes/GML001-Q1851/Brad/seqFISH/svz.h5ad
    OUTPUT: * /Volumes/GML001-Q1851/Brad/seqFISH/svz_gridded.h5ad
"""

################################################################################
                     # Environment Setup #
################################################################################
# TODO change this to your directory
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/stlearn_manuscript/mainfigCCI_newCCISupps/'
import os
os.chdir(work_dir)

import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

import stlearn as st

data_dir = '/Volumes/GML001-Q1851/Brad/seqFISH/'
out_dir = data_dir

################################################################################
                     # Loading the data #
################################################################################
data = sc.read_h5ad(data_dir+'svz.h5ad')

################################################################################
                     # Gridding the data #
################################################################################
# Just example to show comparison with original single cell data #
grid_data = st.tl.cci.grid(data, n_row=30, n_col=30, use_label='cell_type')
sc.pl.spatial(grid_data, color='Gas6', vmax=5)
sc.pl.spatial(grid_data, color='n_cells')

# Actual gridding where reduce resolution #
grid_data = st.tl.cci.grid(data, n_row=15, n_col=10, use_label='cell_type')
print(sum(grid_data.obs['n_cells']), len(data))
sc.pl.spatial(grid_data, color='Gas6', vmax=5, size=2)
sc.pl.spatial(grid_data, color='n_cells', size=2)
sc.pl.spatial(data, color='grid', legend_loc=None)

max_indices = np.apply_along_axis(np.argmax, 1, grid_data.uns['cell_type'].values)
cell_set = np.unique(grid_data.uns['cell_type'].columns.values)
grid_data.obs['cell_type'] = [cell_set[index] for index in max_indices]
sc.pl.spatial(grid_data, color='cell_type')

color_map = {}
for i, ct in enumerate(data.obs['cell_type'].cat.categories):
    color_map[ct] = data.uns['cell_type_colors'][i]

grid_data.uns['cell_type_colors'] = [color_map[ct]
                            for ct in grid_data.obs['cell_type'].cat.categories]

grid_data.obsm["deconvolution"] = grid_data.uns['cell_type']
st.pl.deconvolution_plot(grid_data, show_donut=False, spot_size=55,
                         colors=grid_data.uns['cell_type_colors'])

################################################################################
                     # Running LR analysis #
################################################################################
# Loading the LR databases available within stlearn (from NATMI)
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='mouse')
print(len(lrs))

# Running the analysis #
st.tl.cci.run(grid_data, lrs,
                  min_spots = 3, #Filter out any LR pairs with no scores for less than min_spots
                  distance=0, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=10000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=2, # Number of CPUs for parallel. If None, detects & use all available.
                  )

lr_info = grid_data.uns['lr_summary'] # A dataframe detailing the LR pairs ranked by number of significant spots.
print('\n', lr_info)

st.pl.lr_result_plot(grid_data, 'Gas6_Axl', 'lr_sig_scores', size=300)
plt.show()

# Looks very similar to the original results #

################################################################################
                     # Running the CCI analysis #
################################################################################
# Running the counting of co-occurence of cell types and LR expression #
st.tl.cci.run_cci(grid_data, 'cell_type', # Spot cell information either in data.obs or data.uns
                  min_spots=1, # Minimum number of spots for LR to be tested.
                  spot_mixtures=True, # If True will use the label transfer scores,
                                      # so spots can have multiple cell types if score>cell_prop_cutoff
                  cell_prop_cutoff=0,
                  sig_spots=True, # Only consider neighbourhoods of spots which had significant LR scores.
                  n_perms=1000, # Permutations of cell information to get background, recommend ~1000
                 )

st.pl.lr_chord_plot(grid_data, 'cell_type', 'Gas6_Axl')

### Saving the results #####
#grid_data.write_h5ad(out_dir+'svz_gridded.h5ad', compression='gzip')
grid_data.write_h5ad(out_dir+'svz_gridded2.h5ad', compression='gzip')



"""
Generates a uniform background of randomly selected pairs.

            INPUT: * /Volumes/GML001-Q1851/Brad/skin_analysis/skin_ppd.h5ad
            OUTPUT: * /Volumes/GML001-Q1851/Brad/skin_analysis/naiive_bg.pkl
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

import scripts.X3_skin_analysis.helpers as fhs
from stlearn.tools.microenv.cci.base import calc_neighbours, calc_distance
import scripts.utils.load_data.simple_pickle as spl

data_dir = '/Volumes/GML001-Q1851/Brad/skin_analysis/'
out_dir = data_dir

################################################################################
                            # Load the data  #
################################################################################
data = sc.read_h5ad(data_dir+'skin_ppd.h5ad')

################################################################################
                    # Generating the naiive background  #
################################################################################
# Making sure does not contain any LR pairs #
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit', 'connectomeDB2020_put'])
lr_genes = np.unique([lr_.split('_') for lr_ in lrs])
remaining = [gene for gene in data.var_names if gene not in lr_genes]

data = data[:,remaining]

# Need to calculate the neighbourhoods, same as calculated for other runs. #
distance = calc_distance(data, None)
neighbours = calc_neighbours(data, distance, verbose=True)
data.obsm['spot_neighbours'] = pd.DataFrame([','.join(x.astype(str))
                                              for x in neighbours],
                                             index=data.obs_names,
                                             columns=['neighbour_indices'])

# Generating the background #
background, rand_pairs = fhs.get_stats_naiive_bg(data, 10000, get_stats=False)

################################################################################
                        # Saving the results  #
################################################################################
naiive = {'background': background, 'rand_pairs': rand_pairs}
spl.saveAsPickle(out_dir+'naiive_bg.pkl', naiive)


"""
 Runs code for performing the hotspoting approach to the LR analysis.

INPUT:

    * GML001-Q1851\Jon\Human_Breast_Cancer_Block_A_Section_1

"""

# TODO: NOTE this takes a while to run (~8.5hrs), was run with 23 CPUs on tinaroo.

################################################################################
                # Environment setup & load test data #
################################################################################
#TODO: NOTE must be run in the base directory: stlearn_reproduce_results/
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/SpatialTranscriptomics/' \
           'stlearn_reproduce_results/' # TODO update this with your path
gml_rdm_path = '/Volumes/GML001-Q1851/'
out_dir = gml_rdm_path+'Brad/'

# Assuming running from tinaroo home directory #
import os, sys
os.chdir(work_dir)
import scripts.utils.helpers as hs

import numba
import numpy as np
from timeit import time
from sklearn.decomposition import PCA
import stlearn as st
from stlearn.tools.microenv.cci.base import calc_neighbours, get_lrs_scores, \
                                                                   calc_distance

data_set = 'breast'
out_dir = out_dir+data_set+'/'
jon_dir = gml_rdm_path+'Jon/'
bdata = hs.load_breast( jon_dir )

################################################################################
          # Preselecting LR pairs that show most variation in data #
################################################################################
# TODO: NOTE this ensures faster run-time.

#### Dry run to get candidate LRs ####
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'])
#lrs = np.random.choice(lrs, 50, replace=False) # for testing

min_spots = 20
distance = None
mode = '' if type(distance)==type(None) else 'within-spot_'

### Going to try filtering pairs further by
### doing a PCA on the LR pairs & getting top X pairs across pairs
### from the top PCs.
distance = calc_distance(bdata, distance)
neighbours = calc_neighbours(bdata, distance, verbose=True)
het_vals = np.array([1] * len(bdata))

lr_scores, lrs = get_lrs_scores(bdata, lrs, neighbours, het_vals, 0)

n_ = 20
pca = PCA(n_components=n_)
pca.fit(lr_scores)
pc_sums = np.array([sum(pca.explained_variance_ratio_[0:i])
                for i in range(n_)])
top_pcs = pc_sums < .6
top_pcs = top_pcs if np.any(top_pcs) else [True]*n_
pc_weights = np.abs(pca.components_[top_pcs,:])
totals = pc_weights.sum(axis=1)
pc_weights = np.apply_along_axis(np.divide, 0, pc_weights, totals)
# Getting lrs with top weighting #
cut = .9
lrs_select = []
for i in range(pc_weights.shape[0]):
        x = pc_weights[i,:]
        order = np.argsort(-x)
        for j in range(1, pc_weights.shape[1]):
                weight_sum = sum(x[order[0:j]])
                if weight_sum > .9:
                        lrs_select.extend(lrs[order[0:j]])
                        break
lrs_select = np.unique(lrs_select)
# selected only 30% of the lrs which explained ~60% of variation in data.
lrs = lrs_select

################################################################################
                        # Running stLearn CCI #
################################################################################
n_pairs = 10000
n_cpus = 23
numba.set_num_threads(n_cpus) # Setting no. of threads to run in parallel, set to same as local for timing.
start = time.time()
st.tl.cci.run(bdata, lrs,
                  min_spots = min_spots, #Filter out any LR pairs with no scores for less than 5 spots
                  distance=distance, #40, #distance=0 for within-spot mode
                  n_pairs=n_pairs, #Number of random pairs to generate
                  adj_method='fdr_bh', #MHT correction method
                  min_expr=0, #min expression to be considered as scored.
                  pval_adj_cutoff=.05, save_bg=False,
)
end = time.time()
diff = end-start
print(f'{diff} secs')
print(f'{diff/60} mins')

# Original version #
# bdata.write(out_dir+f'{data_set}_{mode}bgPerLR_lrsubset_spotsubset_noBgs.h5ad',
#                                                              compression='gzip')

# For reproducibility version #
bdata.write(out_dir+f'{data_set}_{mode}'
                    f'bgPerLR_lrsubset_spotsubset_noBgs_reproduced.h5ad',
                                                             compression='gzip')









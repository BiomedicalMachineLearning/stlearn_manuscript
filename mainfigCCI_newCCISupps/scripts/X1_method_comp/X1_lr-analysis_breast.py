""" Runs code for performing the hotspoting approach to the LR analysis.

INPUT:

    * GML001-Q1851\Jon\Human_Breast_Cancer_Block_A_Section_1
    * GML001-Q1851\Jon\C1_skin_raw

"""

################################################################################
                # Environment setup & load test data #
################################################################################
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/'
gml_rdm_path = '/Volumes/GML001-Q1851/'
out_plots = 'plots/bg_eval/'
out_dir = 'data/bg_eval/'

# Assuming running from tinaroo home directory #
import scripts.helpers as hs

import os, sys
os.chdir(work_dir)

import numba
import numpy as np
import stlearn as st
from timeit import time

data_sets = ['breast', 'skin']
data_set = data_sets[0]
out_dir = out_dir+data_set+'/'
jon_dir = gml_rdm_path+'Jon/'
if data_set == 'breast':
    bdata = hs.load_breast( jon_dir )
else:
    bdata = hs.load_skin( jon_dir )

#### Dry run to get candidate LRs ####
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit', 'connectomeDB2020_put'])
lrs = np.random.choice(lrs, 1000, replace=False)

n_pairs = 10000
start = time.time()
st.tl.cci.run(bdata, lrs,
                  use_label = None, #Need to add the label transfer results to object first, above code puts into 'label_transfer'
                  use_het = 'cell_het', #Slot for cell het. results in adata.obsm, only if use_label
                  min_spots = 1, #Filter out any LR pairs with no scores for less than 5 spots
                  distance=None, #40, #distance=0 for within-spot mode
                  n_pairs=n_pairs, #Number of random pairs to generate
                  adj_method='fdr_bh', #MHT correction method
                  lr_mid_dist = 50, #Controls how LR pairs grouped when creating bg distribs, higher number results in less groups
                                    #Recommended to re-run a few times with different values to ensure results robust to this parameter.
                  min_expr=0, #min expression to be considered as scored.
                  pval_adj_cutoff=.05,
)
end = time.time()
diff = end-start
print(f'{diff} secs')
print(f'{diff/60} mins')

bdata.write(out_dir+f'{data_set}_bgPerLR_spotsubset.h5ad', compression='gzip')

# NOTE this was far too slow, even when subsetting, so running on Tinaroo instead.


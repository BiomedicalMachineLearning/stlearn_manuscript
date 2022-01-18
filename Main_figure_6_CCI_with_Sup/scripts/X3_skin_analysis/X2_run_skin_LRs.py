"""
Loads in the preprocessed skin data & performs the LR analysis on this using
    literature reviewed LRs.
            NOTE: run on HPC, should be reasonably fast locally though.

            INPUT:  * /Volumes/GML001-Q1851/Brad/skin_analysis/skin_ppd.h5ad
            OUTPUT: * /Volumes/GML001-Q1851/Brad/skin_analysis/skin_lrs.h5ad
"""

################################################################################
                        # Environment setup  #
################################################################################
# TODO: NOTE must be run in folder with README entitled StLearn Reproduce
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/'
import os
os.chdir(work_dir)

import scanpy as sc
import stlearn as st
from timeit import time

data_dir = '/Volumes/GML001-Q1851/Brad/skin_analysis/'
out_dir = data_dir

################################################################################
                        # Load the data  #
################################################################################
data = sc.read_h5ad(data_dir+'skin_ppd.h5ad')

################################################################################
                   # Running the LR analysis #
################################################################################
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'])

n_pairs = 10000 #100000
n_cpus = 3
start = time.time()
st.tl.cci.run(data, lrs,
                  min_spots = 2, #Filter LR pairs with no scores less than xspots
                  distance=None, #40, #distance=0 for within-spot mode
                  n_pairs=n_pairs, #Number of random pairs to generate
                  adj_method='fdr_bh', #MHT correction method
                  min_expr=0, # Min expression to be considered as scored.
                  pval_adj_cutoff=.05, save_bg=True, n_cpus=n_cpus
)
end = time.time()
diff = end-start
print(f'{diff} secs')
print(f'{diff/60} mins')

#data.write_h5ad(out_dir+'skin_lrs.h5ad', compression='gzip')







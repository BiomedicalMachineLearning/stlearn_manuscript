"""
Running the stlearn LR-CCI analyses on the seqFISH data.

         INPUT: * /Volumes/GML001-Q1851/Brad/seqFISH/svz.h5ad
         OUTPUT: * /Volumes/GML001-Q1851/Brad/seqFISH/svz_LR-CCI.h5ad
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
import matplotlib.pyplot as plt

import stlearn as st

data_dir = '/Volumes/GML001-Q1851/Brad/seqFISH/'
out_dir = data_dir

################################################################################
                     # Loading the data #
################################################################################
data = sc.read_h5ad(data_dir+'svz.h5ad')

################################################################################
                     # Running the LR analysis #
################################################################################
# Loading the LR databases available within stlearn (from NATMI)
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='mouse')
print(len(lrs))

# Running the analysis #
st.tl.cci.run(data, lrs,
                  min_spots = 3, #Filter out any LR pairs with no scores for less than min_spots
                  distance=145, # None defaults to spot+immediate neighbours; distance=0 for within-spot mode
                  n_pairs=10000, # Number of random pairs to generate; low as example, recommend ~10,000
                  n_cpus=2, # Number of CPUs for parallel. If None, detects & use all available.
                  )

lr_info = data.uns['lr_summary'] # A dataframe detailing the LR pairs ranked by number of significant spots.
print('\n', lr_info)

################################################################################
                     # Running the CCI analysis #
################################################################################
# Running the counting of co-occurence of cell types and LR expression #
st.tl.cci.run_cci(data, 'cell_type', # Spot cell information either in data.obs or data.uns
                  min_spots=3, # Minimum number of spots for LR to be tested.
                  spot_mixtures=False, # If True will use the label transfer scores,
                                      # so spots can have multiple cell types if score>cell_prop_cutoff
                  sig_spots=True, # Only consider neighbourhoods of spots which had significant LR scores.
                  n_perms=1000 # Permutations of cell information to get background, recommend ~1000
                 )

################################################################################
               # Saving the results for visualisation later #
################################################################################
data.write_h5ad(out_dir+'svz_LR-CCI.h5ad', compression='gzip')

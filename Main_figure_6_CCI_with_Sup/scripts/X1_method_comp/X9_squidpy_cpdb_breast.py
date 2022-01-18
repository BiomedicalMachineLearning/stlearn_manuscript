"""
Trying the squidpy implimentation of cpdb, since from the API
            (https://squidpy.readthedocs.io/en/latest/api/squidpy.gr.ligrec.html#squidpy.gr.ligrec)
            can clearly update the database. Following tutorial here:
            https://squidpy.readthedocs.io/en/latest/auto_examples/graph/compute_ligrec.html

            INPUT: * /Volumes/GML001-Q1851/Brad/breast_LR&CCIResults.h5ad
            OUTPUT: * data/breast/squidpy/*
"""

################################################################################
                         # Environment setup #
################################################################################
#TODO: NOTE must be run in the base directory: stlearn_reproduce_results/
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/SpatialTranscriptomics/' \
           'stlearn_reproduce_results/' # TODO update this with your path
gml_rdm_path = '/Volumes/GML001-Q1851/'

import os, sys
os.chdir(work_dir)

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import stlearn as st
import squidpy as sq

data_dir = '/Volumes/GML001-Q1851/Brad/'
out_dir = 'data/breast/squidpy/'

################################################################################
                        # Loading the data #
################################################################################
data = sc.read_h5ad(data_dir+'breast_LR&CCIResults.h5ad')

# Getting LRs to test #
lrs = data.uns['lr_summary'].index.values.astype(str)

# Formatting for squidpy input, as per API:
# https://squidpy.readthedocs.io/en/latest/api/squidpy.gr.ligrec.html#squidpy.gr.ligrec
source = [lr.split('_')[0] for lr in lrs]
target = [lr.split('_')[1] for lr in lrs]
int_dict = {'source': source, 'target': target}

################################################################################
                        # Running squidpy #
################################################################################
result = sq.gr.ligrec(data, n_perms=1000, cluster_key="cell_type",
                      copy=True, use_raw=False, interactions=int_dict,
                      show_progress_bar = True,
                      corr_method='fdr_bh', corr_axis='clusters',
                    )
pvals = pd.DataFrame(result['pvalues'])

################################################################################
        # Creating a CCI summary equivalent to other methods #
################################################################################
lrs_ = []
n_cci_sig = []
ccis = []
ps = []
cci_names = np.array([col[0]+'--'+col[1] for col in pvals.columns])
for i, row in enumerate(pvals.index):
    lr_ = '_'.join(list(row))
    lrs_.append(lr_)

    # Getting sig CCIs for this lr #
    lr_pvals = np.array(pvals.values[i,:])
    sig_bool = lr_pvals < .05
    lr_ccis = cci_names[sig_bool]
    lr_sig_pvals = lr_pvals[sig_bool]
    lr_sig_pvals[lr_sig_pvals==0] = 1/1000 #can't be 0!
    lr_sig_logpvals = -np.log10(lr_sig_pvals)

    # Adding to summary output #
    n_cci_sig.append( len(lr_sig_pvals) )
    ccis.append( ','.join(lr_ccis) )
    ps.append( ','.join(lr_sig_logpvals.astype(str)) )

method='squidpy'
features = [f'n_cci_sig_{method}', f'ccis_{method}', f'ps_{method}']
cci_summary = pd.DataFrame([n_cci_sig, ccis, ps]).transpose()
cci_summary.index = lrs_
cci_summary.columns = features

cci_summary.to_csv(out_dir+'squidpy_cci_summary.txt', sep='\t')


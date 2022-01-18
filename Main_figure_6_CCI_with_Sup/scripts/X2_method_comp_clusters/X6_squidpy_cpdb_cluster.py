"""
Trying the squidpy implimentation of cpdb, since from the API
Following tutorial here:
https://squidpy.readthedocs.io/en/latest/auto_examples/graph/compute_ligrec.html

            INPUT: * /Volumes/GML001-Q1851/Brad/breast_ClusterCCIResults.h5ad
            OUTPUT: * data/breast/cluster/
"""

################################################################################
                         # Environment setup #
################################################################################
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/'

import os
os.chdir(work_dir)

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

data_dir = '/Volumes/GML001-Q1851/Brad/'
out_dir = 'data/breast/cluster/'

################################################################################
                        # Loading the data #
################################################################################
data = sc.read_h5ad(data_dir+'breast_ClusterCCIResults.h5ad')

################################################################################
                        # Running squidpy #
################################################################################
result = sq.gr.ligrec(data, n_perms=1000, cluster_key="cell_type",
                      copy=True, use_raw=False,
                      show_progress_bar = True,
                      corr_method='fdr_bh', corr_axis='clusters',
                    )
pvals = pd.DataFrame(result['pvalues'])

################################################################################
     # Creating a CCI interaction matrix equivalent to other methods #
################################################################################
cci_names = np.array([col[0]+'--'+col[1] for col in pvals.columns])
cell_type_set = np.unique([col[0] for col in pvals.columns])
int_matrix = np.zeros((len(cell_type_set), len(cell_type_set)))
for i, row in enumerate(pvals.index):
    lr_ = '_'.join(list(row))

    # Getting sig CCIs for this lr #
    lr_pvals = np.array(pvals.values[i,:])
    sig_bool = lr_pvals < .05
    lr_ccis = cci_names[sig_bool]
    for j, cci in enumerate(lr_ccis):
        c1, c2 = cci.split('--')
        row = np.where(cell_type_set==c1)[0][0]
        col = np.where(cell_type_set==c2)[0][0]
        int_matrix[row,col] += 1

int_df = pd.DataFrame(int_matrix, index=cell_type_set, columns=cell_type_set)

int_df.to_csv(out_dir+'squidpy_ints.txt', sep='\t', header=True)


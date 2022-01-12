"""
Prepares the breast data to be run with cellphonedb, then runs for
                the breast cluster annots.

        Refer here for Laura reference on cellphonedb:
        https://github.com/ellefeg/RNASeq/blob/master/scRNASeq_101/Seurat2CellphoneDB.md
        Also here for tutorial:
        https://github.com/Teichlab/cellphonedb
        Apparently inputting pre-normalised data:
        https://github.com/Teichlab/cellphonedb/issues/12

        INPUT: * /Volumes/GML001-Q1851/Brad/breast_ClusterCCIResults.h5ad
        OUTPUT: * data/breast/cluster/cellphonedb/*
                * data/breast/cluster/cellphonedb/db/*
                                                  -> Running native database
"""

################################################################################
                        # Environment setup #
################################################################################
# TODO: NOTE must be run in folder with README entitled StLearn Reproduce
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

out_dir = 'data/breast/cluster/cellphonedb/'
out_dir2 = 'data/breast/cluster/cellphonedb/db/' #Running with native db

################################################################################
 # Loading data, adding label transfer, & calculate cell type interactions ! #
################################################################################
# Loading the data #
data = sc.read_h5ad(data_dir+f'breast_ClusterCCIResults.h5ad')

labels = data.obs['cell_type'].values.astype(str)

################################################################################
        # Saving the cell type information as per cpdb example data #
################################################################################
cell_data = pd.DataFrame(labels, columns=['cell_type'])
cell_data['Cell'] = data.obs_names
cell_data = cell_data.iloc[:,[1,0]]
cell_data.values[:,1] = [f'_{label}_' for label in cell_data.values[:,1]]
print(cell_data.shape)
cell_data.to_csv(out_dir+'breast_spot_labels.txt', sep='\t', index=False)

################################################################################
                        # Running CellPhoneDB #
################################################################################
os.system(f'cellphonedb method statistical_analysis '
          f'{out_dir}breast_spot_labels.txt '
          f'{data_dir}breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad '
          f'--output-path {out_dir2} --counts-data gene_name')

################################################################################
            # Summarising the CellPhoneDB CCI results #
################################################################################
cpdb_results = pd.read_csv(out_dir2+'pvalues.txt', sep='\t') # 1099 LR pairs.

valid_cci = ['|' in col for col in cpdb_results.columns]
res_sub = cpdb_results.loc[:,valid_cci]

cell_names = np.unique([cci.split('|')[0] for cci in res_sub.columns])
cci_names = np.array([cci.replace('|', '--') for cci in res_sub.columns])

int_matrix = np.zeros((len(cell_names), len(cell_names)))
for i in range(res_sub.shape[0]):
    pvals = res_sub.values[i,:]
    sig = pvals < .05
    cci_sig = cci_names[sig]
    for j, cci in enumerate(cci_sig):
        c1, c2 = cci.split('--')
        row = np.where(cell_names==c1)[0][0]
        col = np.where(cell_names==c2)[0][0]
        int_matrix[row,col] += 1

cell_names2 = [name.replace('_', '') for name in cell_names]
int_df = pd.DataFrame(int_matrix, index=cell_names2, columns=cell_names2)

int_df.to_csv(out_dir+'../cellphonedb_ints.txt', sep='\t')


"""
Prepares the breast data to be run with
                    cellphonedb, then runs for the breast spatial simulated data.

            INPUT: * data/sim_data/spatialsim_v2.h5ad
            OUTPUT: * data/sim_data/methods_out/cellphonedb/*
                    * data/sim_data/methods_out/cellphonedb_ints.txt
"""

################################################################################
                        # Environment setup #
################################################################################
# TODO update these paths
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/'
stlearn_path = '/Users/uqbbalde/Desktop/Uni_Studies/myPython/' \
               'stlearn_latest/stLearn/'
data_dir = 'data/sim_data/s'

import os, sys
os.chdir(work_dir)
sys.path.append(stlearn_path)

import numpy as np
import pandas as pd
import scanpy as sc

out_dir = 'data/sim_data/methods_out/cellphonedb/'

################################################################################
 # Loading data, adding label transfer, & calculate cell type interactions ! #
################################################################################
# Loading the data #
data = sc.read_h5ad(data_dir+f'spatialsim_v2.h5ad')

labels = data.obs['cell_type'].values.astype(str)

################################################################################
        # Saving the cell type information as per example data #
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
run_str = f'cellphonedb method statistical-analysis '+\
          f'{out_dir}breast_spot_labels.txt '+\
          f'{data_dir}spatialsim_v2.h5ad '+\
          f'--output-path {out_dir} --counts-data gene_name'
print(run_str)
os.system(run_str)
# NOTE had to add the cellphonedb code to get this to run !!!!

################################################################################
            # Summarising the CellPhoneDB CCI results #
################################################################################
cpd_int_results = pd.read_csv(out_dir+'interaction_result.txt', sep='\t')
cpd_pvals = pd.read_csv(out_dir+'result_percent.txt', sep='\t')

index_match = [np.where(cpd_int_results.values[:,0]==index)[0][0]
               for index in cpd_pvals.values[:,0]]
cpd_int_results = cpd_int_results.iloc[index_match]
cpd_int_results.index = cpd_pvals.index

cpdb_results = pd.concat([cpd_int_results, cpd_pvals], axis=1) #939 LR pairs

#cpdb_results = pd.read_csv(out_dir+'pvalues.txt', sep='\t') # 1099 LR pairs.
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
cell_names2 = [name if name!='luminalar' else 'luminal_ar'
               for name in cell_names2]
int_df = pd.DataFrame(int_matrix, index=cell_names2, columns=cell_names2)

int_df.to_csv(out_dir+'../cellphonedb_ints.txt', sep='\t')


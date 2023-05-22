"""
Runs the NATMI pipeline on the simulated data.

        INPUT: * data/sim_data/spatialsim_v2.h5ad
        OUTPUT: * data/sim_data/methods_out/natmi/*
                * data/sim_data/methods_out/natmi_ints.txt
"""

################################################################################
                        # Environment setup #
################################################################################
# TODO update to your work dir
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

import matplotlib.pyplot as plt

import scripts.utils.visualisation.quick_plots as qpl
qpl.setUp()

out_dir = 'data/sim_data/methods_out/natmi/'
out_dir2 = 'data/sim_data/methods_out/'

################################################################################
        # Loading data, & saving in format for input with NATMI #
################################################################################
# Loading the data #
data = sc.read_h5ad(data_dir+f'spatialsim_v2.h5ad')

# Writing out in compatible input according to docs here:
#                                           https://github.com/forrest-lab/NATMI
expr = data.to_df().transpose()
expr.columns = [f'bc{i}' for i in expr.columns]
expr.to_csv(out_dir+'breast_expr.csv',index=True,header=True)
ct = data.obs['cell_type']
ct = pd.DataFrame([f'_{i}_' for i in ct], index=expr.columns, columns=['cell_type'])
ct = ct['cell_type']
ct.to_csv(out_dir+'metadata.csv',index=True,header=True)

################################################################################
                            # Running NATMI #
################################################################################
# TODO need to git clone NATMI code from here:   https://github.com/asrhou/NATMI
""" Below we:
        1. cd to the git clone NATMI folder.
        2. activate conda environment
        3. run NATMI on the simulated data.
                   
Run:

cd NATMI/
source activate STI
python ExtractEdges.py --interDB lrdbs/connectomeDB2020_lit.tsv --emFile path/data/sim_data/methods_out/natmi/breast_expr.csv --annFile path/data/sim_data/methods_out/natmi/metadata.csv --coreNum 2 --out path/data/sim_data/methods_out/natmi/
"""

################################################################################
               # Geting the CCIs and N-CCIs per LR #
################################################################################
natmi_res = pd.read_csv(out_dir+'Edges_connectomeDB2020_lit.csv')
ls = natmi_res.loc[:,'Ligand symbol'].values.astype(str)
rs = natmi_res.loc[:,'Receptor symbol'].values.astype(str)
natmi_lrs = np.array([ls[i]+'_'+rs[i] for i in range(len(ls))])

# Getting the ccis #
send = natmi_res.loc[:,'Sending cluster'].values.astype(str)
target = natmi_res.loc[:,'Target cluster'].values.astype(str)
natmi_ccis = np.array(['--'.join([send[i], target[i]])
                       for i in range(len(send))])

# Plotting the expression weights filtering strategy #
qpl.distrib(natmi_res.loc[:,'Edge average expression weight'].values, bins=100,
            log=True)
fig, ax = qpl.distrib(natmi_res.loc[:,
                          'Edge average expression derived specificity'].values,
            bins=100, log=False, show=False)
ax.vlines(.09, 0, 8000)
plt.show()

# Retrieving the LR cci counts #
col = 'Edge average expression derived specificity'
lr_scores = natmi_res.loc[:,col].values
cut = .09 # Based on the above...
unique_lrs = np.unique(natmi_lrs)
cell_type_set = np.unique(list(send)+list(target))
cci_ints = np.zeros((len(cell_type_set), len(cell_type_set)))
for i, lr in enumerate(unique_lrs):
    lr_bool = np.logical_and(natmi_lrs==lr, lr_scores>cut)
    lr_results = natmi_res.loc[lr_bool,:]
    for j in range(lr_results.shape[0]):
        c1 = lr_results.iloc[j,:]['Sending cluster']
        c2 = lr_results.iloc[j,:]['Target cluster']
        row = np.where(cell_type_set==c1)[0][0]
        col = np.where(cell_type_set==c2)[0][0]
        cci_ints[row, col] += 1

cell_type_set_ = [label.strip('_') for label in cell_type_set]
int_df = pd.DataFrame(cci_ints, index=cell_type_set_, columns=cell_type_set_)

# Saving the output!!!! #
int_df.to_csv(out_dir+'../natmi_ints.txt', sep='\t', header=True)



"""
 Runs the NATMI pipeline.

                INPUT: * /Volumes/GML001-Q1851/Brad/breast_LR&CCIResults.h5ad
                OUTPUT: * data/method_comp/breast/natmi/
"""

################################################################################
                        # Environment setup #
################################################################################
#TODO: NOTE must be run in the base directory: stlearn_reproduce_results/
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/SpatialTranscriptomics/' \
           'stlearn_reproduce_results/' # TODO update this with your path
gml_rdm_path = '/Volumes/GML001-Q1851/'
data_dir = gml_rdm_path+'Brad/'

import os, sys
os.chdir(work_dir)

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import scripts.utils.visualisation.quick_plots as qpl
qpl.setUp()

data_set = 'breast'
out_dir = 'data/breast/natmi/'

################################################################################
        # Loading data, & saving in format for input with NATMI #
################################################################################
# Loading the data #
data = sc.read_h5ad(data_dir+f'breast_LR&CCIResults.h5ad')

# Writing out in compatible input according to docs here:
#                                           https://github.com/forrest-lab/NATMI
data.to_df().transpose().to_csv(out_dir+'breast_expr.csv',index=True,header=True)
data.obs['cell_type'].to_csv(out_dir+'metadata.csv',index=True,header=True)

# The st_lrs we are interested in #
st_lrs = list(data.uns['per_lr_cci_cell_type'].keys())

################################################################################
                            # Running NATMI #
################################################################################
""" NATMI can only be run from the directory where it is git cloned from!
    git clone this repo, then cd on the commandline to that directory: 
        https://github.com/asrhou/NATMI 

Run:

# TODO need to update this with the path where this github installed!
# i.e. to where slearn_reproduce_results/ is.
python ExtractEdges.py --interDB lrdbs/connectomeDB2020_lit.tsv --emFile path/data/breast/natmi/breast_expr.csv --annFile path/data/method_comp/breast/natmi/metadata.csv --coreNum 2 --out path/data/breast/natmi/
"""

################################################################################
               # Geting the CCIs and N-CCIs per LR #
################################################################################
natmi_res = pd.read_csv(out_dir+'Edges_connectomeDB2020_lit.csv')
ls = natmi_res.loc[:,'Ligand symbol'].values.astype(str)
rs = natmi_res.loc[:,'Receptor symbol'].values.astype(str)
natmi_lrs = np.array([ls[i]+'_'+rs[i] for i in range(len(ls))])
natmi_bool = [lr in st_lrs for lr in natmi_lrs]
natmi_lrs = natmi_lrs[natmi_bool]
natmi_res = natmi_res.loc[natmi_bool,:]

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
            bins=100, log=True, show=False)
ax.vlines(.015, 0, 8)
plt.show()

# Retrieving the LR cci counts #
col = 'Edge average expression derived specificity'
lr_scores = natmi_res.loc[:,col].values
cut = .015 # Based on the above...
ccis = []
cci_counts = []
cci_scores = []
for i, lr in enumerate(st_lrs):
    lr_bool = np.logical_and(natmi_lrs==lr, lr_scores>cut)
    lr_results = natmi_res.loc[lr_bool,:]
    cci_counts.append( sum(lr_bool) )
    ccis.append( ','.join(natmi_ccis[lr_bool]) )
    cci_scores.append( ','.join(lr_results.loc[:,col].values.astype(str)) )

qpl.distrib(cci_counts, bins=20, log=True, logbase=2)

# Creating the summary dataframe!!!! #
cci_df = pd.DataFrame([cci_counts, ccis, cci_scores], columns=st_lrs,
                      index=['n_cci_sig_NATMI', 'ccis_NATMI', 'ps_NATMI'
                             ]).transpose()
cci_df.to_csv(out_dir+'natmi_cci_summary.txt', sep='\t')



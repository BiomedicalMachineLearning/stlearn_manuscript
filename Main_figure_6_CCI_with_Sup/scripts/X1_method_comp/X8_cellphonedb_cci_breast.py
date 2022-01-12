"""
Prepares the breast data to be run with cellphonedb, then runs.
              Here for tutorial:
              https://github.com/Teichlab/cellphonedb
              Apparently inputting pre-normalised data:
              https://github.com/Teichlab/cellphonedb/issues/12

        INPUT: * /Volumes/GML001-Q1851/Brad/
                                   breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad
               * data/label_transfer_bc.csv
        OUTPUT: * data/method_comp/breast/cellphonedb/*
                * data/method_comp/breast/cellphonedb/db/*
                                    -> Running native database
                * data/method_comp/breast/cellphonedb/st_db/*
                                    -> Running stlearn database
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
import stlearn as st

import scripts.utils.visualisation.quick_plots as qpl
qpl.setUp()

data_set = 'breast'
out_dir = 'data/breast/cellphonedb/'
out_dir4 = 'data/breast/'
data_dir2 = out_dir4

################################################################################
 # Loading data, adding label transfer, & calculate cell type interactions ! #
################################################################################
# Loading the data #
data = sc.read_h5ad(data_dir+
                    f'breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad')

# Adding the label transfer information #
spot_mixtures = pd.read_csv(data_dir2+'label_transfer_bc.csv',
                                                          index_col=0, sep='\t')
labels = spot_mixtures.loc[:,'predicted.id'].values.astype(str)
spot_mixtures = spot_mixtures.drop(['predicted.id','prediction.score.max'],
                                   axis=1)
spot_mixtures.columns = [col.replace('prediction.score.', '')
                         for col in spot_mixtures.columns]
data.obs['cell_type'] = labels

data.uns['cell_type'] = spot_mixtures

################################################################################
        # Saving the cell type information as per example data #
################################################################################
cell_data = pd.DataFrame(labels, columns=['cell_type'])
cell_data['Cell'] = data.obs_names
cell_data = cell_data.iloc[:,[1,0]]
cell_data.to_csv(out_dir+'breast_spot_labels.txt', sep='\t', index=False)

################################################################################
                        # Running CellPhoneDB #
################################################################################
os.system(f'cellphonedb method statistical_analysis '
          f'{out_dir}breast_spot_labels.txt '
          f'{data_dir}breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad '
          f'--output-path {out_dir} --counts-data gene_name')

################################################################################
            # Summarising the CellPhoneDB CCI results #
################################################################################
cpdb_results = pd.read_csv(out_dir+'pvalues.txt', sep='\t') # 1099 LR pairs.

# Subsetting to stlearn lrs #
cci_summary = pd.read_csv(out_dir4+'cci_summary.txt', sep='\t') # 549 LR pairs.
st_lrs = cci_summary.index.values.astype(str)
cpd_lrs = np.array( [cpdb_results.loc[:,'gene_a'].values.astype(str)[i]+\
                       '_'+cpdb_results.loc[:,'gene_b'].astype(str).values[i] \
                       for i in range(cpdb_results.shape[0])] )
cpd_lrs_rev = np.array([cpdb_results.loc[:,'gene_b'].values.astype(str)[i]+\
                       '_'+cpdb_results.loc[:,'gene_a'].astype(str).values[i] \
                       for i in range(cpdb_results.shape[0])])
st_bool = [lr_ in st_lrs for lr_ in cpd_lrs]
st_bool_rev = [lr_ in st_lrs for lr_ in cpd_lrs_rev]
cpdb_sub = cpdb_results.loc[st_bool,:]
print(cpdb_sub.shape[0]) # Only 52 overlapping..
print(cpd_lrs[st_bool])
print(cpd_lrs_rev[st_bool_rev])

# Summarise with the different LRs & get into equivalent format...
cpd_lrs = np.array( [cpdb_results.loc[:,'gene_a'].values.astype(str)[i]+\
                       '_'+cpdb_results.loc[:,'gene_b'].astype(str).values[i] \
                       for i in range(cpdb_results.shape[0])] )
valid_lr_bool = ['nan' not in lr for lr in cpd_lrs]
valid_lrs = cpd_lrs[valid_lr_bool]
valid_cci = ['|' in col for col in cpdb_results.columns]
res_sub = cpdb_results.loc[valid_lr_bool,valid_cci]
res_sub.index = valid_lrs

cci_names = np.array([cci.replace('|', '--') for cci in res_sub.columns])
n_cci_sig = []
ccis = []
ps = []
for i, lr in enumerate(valid_lrs):
    sig_bool = res_sub.values[i,:] < .05
    sig_ps = res_sub.values[i,sig_bool]
    sig_ps[sig_ps==0] = 1/1000
    n_cci_sig.append( len(sig_ps) )
    sig_logps = -np.log10(sig_ps)
    ccis.append( ','.join(cci_names[sig_bool]) )
    ps.append( ','.join(sig_logps.astype(str)) )

summary = pd.DataFrame([n_cci_sig, ccis, ps]).transpose()
summary.index = valid_lrs
method = 'cpdb'
summary.columns = [f'n_cci_sig_{method}', f'ccis_{method}', f'ps_{method}']

summary.to_csv(out_dir+f'{method}_cci_summary.txt', sep='\t')


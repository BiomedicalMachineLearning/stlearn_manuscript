"""
(21/09/21) -> For selected LRs, visualising the
                spatial LR scores along with the cell types as a way to show
                that the stlearn result is correct, where the CellChat
                prediction is incorrect because it lacks the spatial context.

                INPUT:  * /Volumes/GML001-Q1851/Brad/
                               breast_bgPerLR_lrsubset_spotsubset_noBgs.h5ad
                        * data/label_transfer_bc.csv
                OUTPUT: * plots/method_comp/stlearn_breast_spatial_*
"""

################################################################################
                # Environment setup & load test data #
################################################################################
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
import matplotlib.pyplot as plt
import stlearn as st

import stlearn.plotting.cci_plot as cci_plot

import scripts.utils.visualisation.helpers as vhs
import scripts.utils.visualisation.quick_plots as qpl
qpl.setUp()

mode = '' # 'abs_' is discrete spot labels, '' is deconvolution.
data_set = 'breast'
out_plots1 = 'plots/method_comp/'
out_plots = 'plots/method_comp/stlearn_breast_spatial_'

################################################################################
 # Loading data, adding label transfer, & calculate cell type interactions ! #
################################################################################
# Loading the data #
data = sc.read_h5ad(data_dir+'breast_LR&CCIResults.h5ad')

# st.tl.cci_rank.run_cci(data, 'cell_type')

# st.pl.cluster_plot(data, use_label='cell_type', size=10)
# vhs.dealWithPlot(True, True, True, out_plots, 'labels.pdf', 300)

################################################################################
        # Creating the lr result plots for LRs of interest #
# ################################################################################
lrs = ['TGM2_ADGRG1', 'GPC3_IGF1R', 'SPP1_ITGB1']
size = 18
a = .85
for best_lr in lrs:
    show = True
    st.pl.lr_result_plot(data, use_result='lr_scores', use_lr=best_lr,
                         size=size, cell_alpha=a)
    vhs.dealWithPlot(True, show, True, out_plots,
                     f'{best_lr}_lrscores.pdf', 300)

    st.pl.lr_result_plot(data, use_result='lr_sig_scores', use_lr=best_lr,
                         size=size, cell_alpha=a)
    vhs.dealWithPlot(True, show, True, out_plots,
                     f'{best_lr}_lr-sig-scores.pdf', 300)

    # TODO debug below #
    st.pl.lr_plot(data, best_lr, outer_size_prop=1, outer_mode='binary',
                  pt_scale=20, use_label=None, show_image=True, sig_spots=True)
    vhs.dealWithPlot(True, show, True, out_plots,
                     f'{best_lr}_binarySigLR.pdf', 300)

    st.pl.lr_plot(data, best_lr, inner_size_prop=0.1, outer_mode='binary',
                  min_expr=2,
                  pt_scale=10, use_label=None, show_image=True, sig_spots=False)
    vhs.dealWithPlot(True, show, True, out_plots,
                     f'{best_lr}_binaryLR_spatial.png', 300)

###### Just for 'GPC3_IGF1R' adding the arrows ######
best_lr = 'GPC3_IGF1R'
cci_plot.lr_plot(data, best_lr, outer_size_prop=1, outer_mode=None,
              pt_scale=40, use_label='cell_type', show_arrows=True,
              show_image=True, sig_spots=False, sig_cci=True,
                 arrow_head_width=4,
                 arrow_width=1
                 )
vhs.dealWithPlot(True, True, True, out_plots,
                 f'{best_lr}_arrow-cell-types.pdf', 300)

#### Another version with a lower alpha to see structure ######
best_lr = 'GPC3_IGF1R'
cci_plot.lr_plot(data, best_lr, outer_size_prop=1, outer_mode=None,
              pt_scale=40, use_label='cell_type', show_arrows=True,
              show_image=True, sig_spots=False, sig_cci=True,
                 arrow_head_width=4,
                 arrow_width=1, cell_alpha=.8
                 )
vhs.dealWithPlot(True, True, True, out_plots,
                 f'{best_lr}_see-through_arrow-cell-types.pdf', 300)



# Understanding matplotlib arrows #
fig, ax = plt.subplots()
ax.arrow(0, 0, 1, 1, width=.1, head_width=.2,
         linewidth=0
         )
plt.show()



"""
(25/10/21) -> Multi-method comparison of CCI ranking based on no. of significant
                LR pairs. Will focus on pair which is highly ranked in the other
                methods but low-ranked by stLearn as example, and show the
                squidpy colocalisation analysis of the two cell types. Then
                repeat this for every CCI and predicted LR pair in order to plot
                distributions of squidpy spatial AUC to show stLearn enriches
                for cell types spatially when compared to other methods.

              INPUT: * data/method_comp/breast/cci_summary.txt
                     * data/method_comp/breast/singlecellsignalr/
                                                           scsignalr_summary.txt
                     * data/method_comp/breast/natmi/natmi_cci_summary.txt
              OUTPUT: * plots/multi_comp_v2/*
"""

################################################################################
                         # Environment setup #
################################################################################
work_dir = '/Users/uqbbalde/Desktop/Uni_Studies/projects/stLearn/'

import os, sys
os.chdir(work_dir)

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq

import importlib

import scripts.utils.visualisation.helpers as vhs
import scripts.utils.visualisation.quick_plots as qpl
import scripts.utils.multi_comp_helpers as mch

import stlearn.plotting.cci_plot_helpers as cci_hs

data_dir = 'data/method_comp/breast/'
data_dir2 = '/Volumes/GML001-Q1851/Brad/'
data_dir3 = 'data/method_comp/breast/singlecellsignalr/'
data_dir4 = 'data/method_comp/breast/natmi/'
out_plots = 'plots/multi_comp_v2/'

out = '.pdf'

################################################################################
                        # Loading the data #
################################################################################
#### Also loading the LR features to look for expression level bias in CellChat
data = sc.read_h5ad(data_dir2+'breast_LR&CCIResults.h5ad')
lrfeatures = data.uns['lrfeatures']

cci_summary = pd.read_csv(data_dir+'cci_summary.txt', sep='\t')
# Subsetting to cci_rank tested #
cci_summary = cci_summary.loc[cci_summary.loc[:,'cci_tested'], :]

#### Loading the other methods ran ######
signalr_summary = pd.read_csv(data_dir3+'scsignalr_summary.txt', sep='\t')
natmi_summary = pd.read_csv(data_dir4+'natmi_cci_summary.txt',
                            sep='\t', index_col=0)

### merging #####
methods = ['', '_NATMI', '_cellchat', '_scsignalr']
method_colors = ['gold', 'hotpink', 'darkcyan', 'limegreen']
cci_summary = pd.concat([cci_summary, signalr_summary, natmi_summary], axis=1)

################################################################################
         # Plotting the CCIs ranked by no. of LRs significant #
################################################################################
label_set = np.unique(data.obs['cell_type'].values)
cci_set = []
for labeli in label_set:
    for labelj in label_set:
        cci_set.append(f'{labeli}->{labelj}')
cci_set = np.unique(cci_set)

#### Counting the no. of LRs per CCI ######
method_ccis_ranked = []
method_nlrs_ranked = []
method_ccilrs_ranked = []
for i in range(len(methods)):
    ccis_ranked, cci_nlrs_ranked, cci_lrs_ranked = \
                        mch.get_cci_nlr_ranked(cci_set, cci_summary, methods[i])
    method_ccis_ranked.append( ccis_ranked )
    method_nlrs_ranked.append( cci_nlrs_ranked )
    method_ccilrs_ranked.append( cci_lrs_ranked )

# Finding an CCI which has the highest rank by the other methods but
# lowest in stLearn.
perfect = np.array([len(cci_set)] + [0] * (len(methods) - 1))
cand_ccis = mch.get_cands(perfect, cci_set, methods, method_ccis_ranked)

###### Plotting the ranked lists - global view ######
int_cci = cand_ccis[0]
highlight_lrs = [int_cci]
rot = 20
fig, ax = plt.subplots(figsize=(8,3))
for method in methods[::-1]:
    i = np.where(np.array(methods)==method)[0][0]
    ccis_ranked, nlrs_ranked = method_ccis_ranked[i], method_nlrs_ranked[i]
    cci_hs.rank_scatter(ccis_ranked, nlrs_ranked, rot=rot,
                        x_label='CCI Ranked', y_label='n_lrs',
                        color=method_colors[i], ax=ax,
                        highlight_items=highlight_lrs, show_text=True,
                        show=False)
ylims = ax.get_ylim()
ax.set_ylim((-1, 300))
vhs.dealWithPlot(True, True, True,
                 out_plots, f'multicomp_breast_cci_ranked{out}', 300)

################################################################################
     # Getting candidate example CCI & LR
################################################################################
#### Getting the LRs which are called sig for CCI of interest each method ####
int_lrs = []
int_lr_stat = []
for i, method in enumerate(methods):
    # Getting the LRs predicted to interact interesting cell types #
    int_cci_loc = np.where(method_ccis_ranked[i]==int_cci)[0][0]
    method_int_cci_lrs = method_ccilrs_ranked[i][int_cci_loc].split(',')
    remaining_lrs = [lr for lr in cci_summary.index.values
                     if lr not in method_int_cci_lrs]
    method_lrs = np.array( method_int_cci_lrs+remaining_lrs )

    # Getting the sig scores for the lr cci method #
    lr_stats = []
    for j, lr in enumerate(method_int_cci_lrs):
        lr_ccis = np.array(cci_summary.loc[lr,f'ccis{method}'].replace(
                                                         '--', '->').split(','))
        lr_cci_loc = np.where(lr_ccis == int_cci)[0][0]
        lr_cci_stat = float(np.array(cci_summary.loc[lr,
                                         f'ps{method}'].split(','))[lr_cci_loc])
        lr_stats.append( lr_cci_stat )
    lr_stats = np.array( lr_stats + [0]*len(remaining_lrs) )

    # Min-max scaling #
    min_, max_ = min(lr_stats), max(lr_stats)
    stats_scaled = (lr_stats-min_) / (max_ - min_)

    # Ordering and saving #
    order = np.argsort(-stats_scaled)
    method_lrs_ordered = method_lrs[order]
    method_lr_stats_ordered = stats_scaled[order]
    int_lrs.append( method_lrs_ordered )
    int_lr_stat.append( method_lr_stats_ordered )

# Finding an CCI which has the highest rank by the other methods but
# lowest in stLearn.
lr_set = cci_summary.index.values.astype(str)
perfect_low = np.array( [cci_summary.shape[0]] + [0]*(len(methods)-1) )
perfect_high = np.array( [0] + [cci_summary.shape[0]]*(len(methods)-1) )
stlow_cands = mch.get_cands(perfect_low, lr_set, methods, int_lrs)
sthigh_cands = mch.get_cands(perfect_high, lr_set, methods, int_lrs)

###### Plotting the ranked lists of LRs for biggest diff CCI - global view ######
figsize = (8,5)

low_lr = stlow_cands[3]
highlight_lrs = [low_lr]
rot = 10
fig, ax = mch.plot_ranked(rot, methods, int_lrs, int_lr_stat,
                'LRs Ranked', f'{int_cci}', method_colors, highlight_lrs,
                          figsize=figsize)
ylims = ax.get_ylim()
ax.set_ylim((-.1, 1.1))
vhs.dealWithPlot(True, True, True,
               out_plots, f'multicomp_breast_low-lr_{int_cci}_ranked{out}', 300)

##### Now for the high_lr for stlearn ######
high_lr = sthigh_cands[0]
highlight_lrs = [high_lr]
rot = 80
fig, ax = mch.plot_ranked(rot, methods[::-1], int_lrs[::-1], int_lr_stat[::-1],
                'LRs Ranked', f'{int_cci}', method_colors[::-1], highlight_lrs,
                          figsize=figsize)
ylims = ax.get_ylim()
ax.set_ylim((-.1, 1.1))
vhs.dealWithPlot(True, True, True,
               out_plots, f'multicomp_breast_high-lr_{int_cci}_ranked{out}', 300)

####### Trying to plot ranked with histograms on either side of axies ##########
importlib.reload(mch)
highlight_lr = stlow_cands[3]
mch.plot_ranked_wHists(int_lrs, int_lr_stat, methods, highlight_lr,
                       height=8, rot=5)
vhs.dealWithPlot(True, True, True,
         out_plots, f'multicomp_breast_high-lr_{int_cci}_rankedwHist{out}', 300)

highlight_lr = sthigh_cands[0]
mch.plot_ranked_wHists(int_lrs[::-1], int_lr_stat[::-1], methods[::-1],
                       highlight_lr, height=8, rot=20)
vhs.dealWithPlot(True, True, True,
         out_plots, f'multicomp_breast_low-lr_{int_cci}_rankedwHist{out}', 300)

################################################################################
         # Using squidpy, looking at the spatial enrichment between the
         # example sender cell & receiver cell and most significant LR from
         # each method.
################################################################################
##### For the stlow lr #####
l, r = low_lr.split('_')
sender, target = int_cci.split('->')
colors = ['royalblue', 'lime', 'red']
mch.plot_cci(int_cci, low_lr, data, colors)
vhs.dealWithPlot(True, True, True, out_plots,
                 f'breast_stlow-labels_spatial.pdf', 300)

# Performing squidpy enrichment
sq.gr.co_occurrence(data, cluster_key="cci_labels")
sq.pl.co_occurrence(data,
    cluster_key="cci_labels",
    clusters=f'{l}--{sender}',
    figsize=(8, 4),
)
# Getting the max point.. #
label_set = data.obs['cci_labels'].cat.categories.values.astype(str)
i_loc = np.where(label_set == f'{l}--{sender}')[0][0]
j_loc = np.where(label_set == f'{r}--{target}')[0][0]
max_loc = np.argmax(data.uns['cci_labels_co_occurrence']['occ'][i_loc,j_loc,:])
max_val = max(data.uns['cci_labels_co_occurrence']['occ'][2,1,:])
max_dist = data.uns['cci_labels_co_occurrence']['interval'][max_loc]
print(max_val) # 1.3332853
vhs.dealWithPlot(True, True, True, out_plots,
                 f'breast_stlow_co-occurrence.pdf', 300)

###### For the sthigh LR ######
colors = ['royalblue', 'orange', 'aqua', ]
l, r = high_lr.split('_')
mch.plot_cci(int_cci, high_lr, data, colors)
vhs.dealWithPlot(True, True, True, out_plots,
                 f'breast_sthigh-labels_spatial.pdf', 300)
sq.gr.co_occurrence(data, cluster_key="cci_labels")
sq.pl.co_occurrence(data,
    cluster_key="cci_labels",
    clusters=f'{l}--{sender}',
    figsize=(8, 4),
)
# Getting the max point.. #
label_set = data.obs['cci_labels'].cat.categories.values.astype(str)
i_loc = np.where(label_set == f'{l}--{sender}')[0][0]
j_loc = np.where(label_set == f'{r}--{target}')[0][0]
max_loc = np.argmax(data.uns['cci_labels_co_occurrence']['occ'][i_loc,j_loc,:])
max_val = max(data.uns['cci_labels_co_occurrence']['occ'][2,1,:])
max_dist = data.uns['cci_labels_co_occurrence']['interval'][max_loc]
print(max_val) # 3.5628283
vhs.dealWithPlot(True, True, True, out_plots,
                 f'breast_sthigh_co-occurrence.pdf', 300)

################################################################################
    # Calculating the squidpy enrichment values for LR and CCIs as above #
################################################################################
lr_genes = np.unique([lr.split('_') for lr in cci_summary.index])
expr = data.to_df().loc[:,lr_genes]
labels = data.obs['cell_type'].values.astype(str)

### Testing with subset .. ###
n_top = 20
method_spatial_scores = []
n = 0
print( (n*6)/60/60 )
for i, method in enumerate(methods):
    # Getting the top LRs for each method #
    ccis_ranked = method_ccis_ranked[i][0:n_top]

    spatial_scores = []
    for j, cci in enumerate(ccis_ranked):
        sender, target = cci.split('->')

        send_bool = labels == sender
        targ_bool = labels == target

        ccilrs_ranked = method_ccilrs_ranked[i][j].split(',')

        for k, lr in enumerate(ccilrs_ranked[0:n_top]):
            n += 1
            l, r = lr.split('_')
            l_bool = expr.loc[:, l].values.astype(float) > 0
            r_bool = expr.loc[:, r].values.astype(float) > 0

            cci_lr_labels = []
            for z in range(data.shape[0]):
                if send_bool[z] and l_bool[z]:
                    cci_lr_labels.append(f'{l}--{sender}')
                elif targ_bool[z] and r_bool[z]:
                    cci_lr_labels.append(f'{r}--{target}')
                else:
                    cci_lr_labels.append(f'')

            data.obs['cci_labels'] = cci_lr_labels
            data.obs['cci_labels'] = data.obs['cci_labels'].astype('category')

            #### Running squidpy ####
            sq.gr.co_occurrence(data, cluster_key="cci_labels", n_jobs=2)

            #### Retrieving the max spatial enrichment value ###
            label_set = data.obs['cci_labels'].cat.categories.values.astype(str)
            i_loc = np.where(label_set == f'{l}--{sender}')[0][0]
            j_loc = np.where(label_set == f'{r}--{target}')[0][0]
            max_val = max(data.uns['cci_labels_co_occurrence']['occ']
                                                              [i_loc, j_loc, :])
            spatial_scores.append( max_val )

    method_spatial_scores.append( spatial_scores )


##### Plotting the distribution of the spatial enrichment scores #######
fig, ax = plt.subplots(figsize=(8,3))
x = 'Spatial colocalisation of CCI-LRs'
bins = 22
a = .7
log=True
for method in methods[::-1]:
    i = np.where(np.array(methods)==method)[0][0]
    spatial_scores = method_spatial_scores[i]
    qpl.distrib(spatial_scores, ax=ax, fig=fig, color=method_colors[i],
                density=False, label=method, x_label=x, bins=bins,
                                        logbase=2, alpha=a, log=log, show=False)
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles[::-1], labels[::-1])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
vhs.dealWithPlot(True, True, True, out_plots,
                 f'multicomp_breast_spatial-scores-sub_distribs{out}', 300)

# NOTE:
#   1. Full version deployed on tinaroo for the top 10 LR pairs per method with many jobs ASAP.
#   Implimenting this in a new minimalist script that simply save the method_spatial_scores...






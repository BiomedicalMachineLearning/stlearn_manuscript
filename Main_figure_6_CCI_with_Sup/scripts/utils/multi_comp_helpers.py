"""
Helper functions for comparing the results between the methods !!!!
"""

import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
from numba import njit, jit, prange
from numba.typed import List
from collections import defaultdict

from stlearn.tools.microenv.cci.het_helpers import get_neighbourhoods

import stlearn.plotting.cci_plot_helpers as cci_hs
from scipy.spatial.distance import euclidean, canberra

def get_cci_ranked(cci_summary, method):
    """Ranks the LR based on no. of ccis for inputted method.
        Must have column with n_cci_sig{method}
    """
    all_lrs = cci_summary.index.values.astype(str)

    st_ccis = cci_summary.loc[:, f'n_cci_sig{method}'].values
    st_order = np.argsort(-st_ccis)
    st_ccis_ranked = st_ccis[st_order]
    st_lrs_ranked = all_lrs[st_order]
    return st_ccis_ranked, st_lrs_ranked

def get_cci_nlr_ranked(cci_set, cci_summary, method):
    """Gets the CCI ranked by no. of LRs predicted to facilitate communication
    between them.
    """
    cci_counts = np.zeros((len(cci_set)))
    #cci_lrs = [''] * len(cci_set)
    cci_lrs = defaultdict(list)
    cci_lr_stats = defaultdict(list)
    for i, lr in enumerate(cci_summary.index.values):
        lr_cci_str = cci_summary.loc[lr, f'ccis{method}']
        if type(lr_cci_str) != str:
            continue
        lr_ccis = lr_cci_str.replace('--', '->').split(',')
        lr_cci_stats = np.array(cci_summary.loc[lr,
                                        f'ps{method}'].split(',')).astype(float)
        for j, cci in enumerate(lr_ccis):
            cci_loc = np.where(cci_set == cci)[0][0]
            cci_counts[cci_loc] += 1
            # cci_lrs[cci_loc] = cci_lrs[cci_loc] + f',{lr}' \
            #                                if len(cci_lrs[cci_loc])!=0 else lr
            cci_lrs[cci].append( lr )
            cci_lr_stats[cci].append( lr_cci_stats[j] ) #Adding the statistic #

    # Saving & re-ordering the CCI lrs by significance level #
    cci_strs = []
    for cci in cci_set:
        if cci in cci_lrs:
            cci_lrs_ = np.array( cci_lrs[cci] )
            cci_lr_stats_ = np.array( cci_lr_stats[cci] )
            order = np.argsort( -cci_lr_stats_ )
            cci_strs.append( ','.join(cci_lrs_[order]) )
        else: # Has no interaction #
            cci_strs.append( '' )

    cci_strs = np.array(cci_strs)

    order = np.argsort(-cci_counts)
    ccis_ranked = cci_set[order]
    cci_nlrs_ranked = cci_counts[order]
    cci_lrs_ranked = cci_strs[order]

    return ccis_ranked, cci_nlrs_ranked, cci_lrs_ranked

def get_cands(perfect, cci_set, methods, method_ccis_ranked, method_ccis_stats):
    """Gets candidate difference between stLearn and other methods by measuring \
        distance from perfect.
    """
    cci_scores = []
    max_rank = len(method_ccis_stats[0])
    for i, cci in enumerate(cci_set):
        cci_ranks = []
        for j in range(len(methods)):
            if cci not in method_ccis_ranked[j]:
                cci_ranks.append( max_rank )
                continue

            rank = np.where(method_ccis_ranked[j] == cci)[0][0]
            # account for equal rank toward end of list #
            rank = max_rank if method_ccis_stats[j][rank]==0 else rank
            cci_ranks.append( rank )
        score = euclidean(perfect, cci_ranks)
        cci_scores.append(score)
    cand_ccis = cci_set[np.argsort(cci_scores)]

    return cand_ccis

def plot_ranked(rot, methods, int_lrs, int_lr_stat, x_label, y_label,
                method_colors, highlight_lrs, figsize=(10,3)):
    """Plots different rankings on same axis.
    """
    fig, ax = plt.subplots(figsize=figsize)
    for method in methods[::-1]:
        i = np.where(np.array(methods) == method)[0][0]
        lrs_ranked, stat_ranked = int_lrs[i], int_lr_stat[i]
        if type(highlight_lrs)!=type(None):
            highlight_lrs_ = [lr for lr in highlight_lrs if lr in lrs_ranked]
        cci_hs.rank_scatter(lrs_ranked, stat_ranked, rot=rot,
                            x_label= x_label, y_label=y_label,
                            color=method_colors[i], ax=ax,
                            highlight_items=highlight_lrs_, show_text=True,
                            show=False)
    return fig, ax


def plot_ranked_wHists(int_lrs, int_lr_stat, methods, highlight_lr,
                        method_map, color_dict_, height=8, rot=2):
    """ Plotting the ranked data with histograms showing the densities.
    """
    index = ['lrs', 'stats', 'rank', 'method']
    # method_map = {'': 'stLearn', '_NATMI': 'NATMI', '_cellchat': 'CellChat',
    #               '_scsignalr': 'scSignalR'}
    dfs = [
        pd.DataFrame([int_lrs[i], int_lr_stat[i], np.argsort(-int_lr_stat[i]),
                      [method_map[method]] * len(int_lr_stat[i])],
                     index=index).transpose()
        for i, method in enumerate(methods)]
    all_df = pd.concat(dfs)

    # Ensuring same ordering as input methods.. #
    # color_dict_ = {'stLearn': 'gold', 'NATMI': 'hotpink', 'CellChat': 'darkcyan',
    #               'scSignalR': 'limegreen'}
    color_dict = {}
    for method in methods:
        color_dict[method_map[method]] = color_dict_[method_map[method]]

    multivariateGrid('rank', 'stats', 'method', all_df,
                         color_dict=color_dict,
                         show=False, height=height)

    for method in methods:
        method = method_map[method]
        method_loc = all_df.loc[:, 'method'].values.astype(str) == method
        lr_loc = all_df.loc[:, 'lrs'].values.astype(str) == highlight_lr
        loc_bool = np.logical_and(method_loc, lr_loc)
        if not np.any(loc_bool):
            continue

        x, y = all_df.loc[loc_bool, 'rank'], all_df.loc[loc_bool, 'stats']
        plt.scatter([x], [y], alpha=1,
                    c=color_dict[method], s=200, #'red',
                    edgecolors='dimgray',
                    )

        lr_text_fp = {'weight': 'bold', 'size': 12}
        plt.text(x - .2, y, highlight_lr, rotation=rot, fontdict=lr_text_fp
                 )


def multivariateGrid(col_x, col_y, col_k, df, color_dict=None,#k_is_color=False,
                     scatter_alpha=.5, show=True, height=6):
    def colored_scatter(x, y, c=None):
        def scatter(*args, **kwargs):
            args = (x, y)
            if c is not None:
                kwargs['c'] = c
            kwargs['alpha'] = scatter_alpha
            plt.scatter(*args, **kwargs)

        return scatter

    g = sns.JointGrid(
        x=col_x,
        y=col_y,
        data=df, height=height, ylim=(-.23, 1.23), ratio=5,
    )
    color = None
    legends = []
    #for name, df_group in df.groupby(col_k):
    for name in color_dict:
        group_bool = df.loc[:,col_k].values == name
        df_group = df.loc[group_bool,:]
        legends.append(name)
        if color_dict:
            color = color_dict[name]
        g.plot_joint(
            colored_scatter(df_group[col_x], df_group[col_y], color),
        )
        # sns.distplot(
        #     df_group[col_x].values,
        #     ax=g.ax_marg_x,
        #     color=color,
        # )
        sns.distplot(
            df_group[col_y].values,
            ax=g.ax_marg_y,
            color=color, kde=True, hist=False,
            vertical=True
        )
    # Do also global Hist:
    # sns.distplot(
    #     df[col_x].values,
    #     ax=g.ax_marg_x,
    #     color='grey'
    # )
    # sns.distplot(
    #     df[col_y].values.ravel(),
    #     ax=g.ax_marg_y,
    #     color='grey',
    #     vertical=True
    # )
    # plt.legend(legends)
    if show:
        plt.show()

def plot_cci(int_cci, int_lr, data, colors, sig=False):
    """ Plots cci predicted to interact via given lr.
    """
    l, r = int_lr.split('_')
    sender, target = int_cci.split('->')

    # labelling where sender cell expresses ligand, and receiver expresses receptor #
    expr = data.to_df()
    labels = data.obs['cell_type'].values.astype(str)
    l_bool = expr.loc[:, l].values.astype(float) > 0
    r_bool = expr.loc[:, r].values.astype(float) > 0
    both_bool = np.logical_and(l_bool, r_bool)
    send_bool = labels == sender
    targ_bool = labels == target
    if sig: # Subsetting to significant hotspots & neighbours for stLearn #
        neighbours = list(get_neighbourhoods(data)[0])
        lr_index = np.where(data.uns['lr_summary'].index.values==int_lr)[0][0]
        sig_spots = data.obsm['lr_sig_scores'][:, lr_index] > 0
        sig_spot_indices = np.where( sig_spots )[0]
        for index in sig_spot_indices:
            neigh_indices = neighbours[index]
            sig_spots[neigh_indices] = True

        l_bool = np.logical_and(l_bool, sig_spots)
        r_bool = np.logical_and(r_bool, sig_spots)

    cci_lr_labels = []
    for i in range(data.shape[0]):
        # if send_bool[i] and both_bool[i]:
        #     cci_lr_labels.append(f'{l}_{r}--{sender}')
        # elif targ_bool[i] and both_bool[i]:
        #     cci_lr_labels.append(f'{l}_{r}--{target}')
        if send_bool[i] and l_bool[i]:
            cci_lr_labels.append(f'{l}--{sender}')
        # elif send_bool[i]:
        #     cci_lr_labels.append(f'{sender}')
        elif targ_bool[i] and r_bool[i]:
            cci_lr_labels.append(f'{r}--{target}')
        # elif targ_bool[i]:
        #     cci_lr_labels.append(f'{target}')
        else:
            cci_lr_labels.append(f'')

    data.obs['cci_labels'] = cci_lr_labels
    data.obs['cci_labels'] = data.obs['cci_labels'].astype('category')
    data.uns['cci_labels_colors'] = colors

    sc.pl.spatial(data, color='cci_labels', size=1.5, alpha=.8,
                  show=False)

@njit
def get_dists(sources, targets, labels, l_bool, r_bool, spatial):
    """ Getting the distances.
    """
    method_dists = List()
    for k in range(len(sources)):
        source, target = sources[k], targets[k]
        # Getting indices of cell types that have ligand/receptor expression if
        # source/target, respectively.
        source_indices = List()
        target_indices = List()
        for i in range(len(labels)):
            if labels[i] == source and l_bool[i]:
                source_indices.append(i)
            if labels[i] == target and r_bool[i]:
                target_indices.append(i)

        for s in source_indices:
            s_loc = spatial[s, :]
            for t in target_indices:
                t_loc = spatial[t, :]
                dist = np.sqrt(np.sum((s_loc-t_loc)**2))
                method_dists.append( dist )

    return method_dists

################################################################################
                        # Multi-comp Version 3 #
################################################################################
def run_cci_squid(data, cci, labels, ):
    """ Running squidpy colocalisation for given CCI.
    """
    send_bottom, targ_bottom = cci.split('->')
    bottom_labels = labels.copy()
    bottom_labels[np.logical_and(labels!=send_bottom, labels!=targ_bottom)] = ''

    label_key = 'labels'
    data.obs[label_key] = bottom_labels
    data.obs[label_key] = data.obs[label_key].astype('category')

    # Performing squidpy enrichment
    sq.gr.co_occurrence(data, cluster_key=label_key)
    sq.pl.co_occurrence(data,
        cluster_key=label_key,
        clusters=send_bottom,
        figsize=(8, 4),
    )
    # Getting the max point.. #
    label_set = data.obs[label_key].cat.categories.values.astype(str)
    i_loc = np.where(label_set == send_bottom)[0][0]
    j_loc = np.where(label_set == targ_bottom)[0][0]
    max_loc = np.argmax(data.uns[f'{label_key}_co_occurrence']['occ'][i_loc,j_loc,:])
    max_val = max(data.uns[f'{label_key}_co_occurrence']['occ'][i_loc,j_loc,:])
    max_dist = data.uns[f'{label_key}_co_occurrence']['interval'][max_loc]
    print(max_val) # 1.4371388






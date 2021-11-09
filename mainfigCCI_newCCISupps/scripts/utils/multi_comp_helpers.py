"""
Helper functions for comparing the results between the methods !!!!
"""

import numpy as np
import pandas as pd
import seaborn as sns
import scanpy as sc
import matplotlib.pyplot as plt
from numba import njit, jit, prange
from numba.typed import List
from collections import defaultdict

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

def get_cands(perfect, cci_set, methods, method_ccis_ranked):
    """Gets candidate difference between stLearn and other methods by measuring \
        distance from perfect.
    """
    cci_scores = []
    for i, cci in enumerate(cci_set):
        cci_ranks = []
        for j in range(len(methods)):
            cci_ranks.append(np.where(method_ccis_ranked[j] == cci)[0][0])
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
        cci_hs.rank_scatter(lrs_ranked, stat_ranked, rot=rot,
                            x_label= x_label, y_label=y_label,
                            color=method_colors[i], ax=ax,
                            highlight_items=highlight_lrs, show_text=True,
                            show=False)
    return fig, ax


def plot_ranked_wHists(int_lrs, int_lr_stat, methods, highlight_lr,
                       height=8, rot=2):
    """ Plotting the ranked data with histograms showing the densities.
    """
    index = ['lrs', 'stats', 'rank', 'method']
    method_map = {'': 'stLearn', '_NATMI': 'NATMI', '_cellchat': 'CellChat',
                  '_scsignalr': 'scSignalR'}
    dfs = [
        pd.DataFrame([int_lrs[i], int_lr_stat[i], np.argsort(-int_lr_stat[i]),
                      [method_map[method]] * len(int_lr_stat[i])],
                     index=index).transpose()
        for i, method in enumerate(methods)]
    all_df = pd.concat(dfs)

    color_dict = {'stLearn': 'gold', 'NATMI': 'hotpink', 'CellChat': 'darkcyan',
                  'scSignalR': 'limegreen'}
    multivariateGrid('rank', 'stats', 'method', all_df,
                         color_dict=color_dict,
                         show=False, height=height)

    for method in methods:
        method = method_map[method]
        method_loc = all_df.loc[:, 'method'].values.astype(str) == method
        lr_loc = all_df.loc[:, 'lrs'].values.astype(str) == highlight_lr
        loc_bool = np.logical_and(method_loc, lr_loc)

        x, y = all_df.loc[loc_bool, 'rank'], all_df.loc[loc_bool, 'stats']
        plt.scatter([x], [y], alpha=.8, c='red',
                    edgecolors='none',
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
    for name, df_group in df.groupby(col_k):
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

def plot_cci(int_cci, int_lr, data, colors):
    """ Plots cci predicted to interact via given lr.
    """
    l, r = int_lr.split('_')
    sender, target = int_cci.split('->')

    # labelling where sender cell expresses ligand, and receiver expresses receptor #
    expr = data.to_df()
    labels = data.obs['cell_type'].values.astype(str)
    l_bool = expr.loc[:, l].values.astype(float) > 0
    r_bool = expr.loc[:, r].values.astype(float) > 0
    send_bool = labels == sender
    targ_bool = labels == target

    cci_lr_labels = []
    for i in range(data.shape[0]):
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


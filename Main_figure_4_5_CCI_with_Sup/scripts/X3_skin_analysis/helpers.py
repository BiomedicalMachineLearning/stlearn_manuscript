""" Helpers for evaluating background approach.
"""

import numpy as np
from numba.typed import List
from statsmodels.stats.multitest import multipletests

import matplotlib.pyplot as plt

from stlearn.tools.microenv.cci.base import get_lrs_scores
from stlearn.tools.microenv.cci.perm_utils import gen_rand_pairs
import scripts.utils.visualisation.helpers as vhs

def get_stats_naiive_bg(data, n_pairs, get_stats=True):
    """Gets the p-vals & padjs when using a naiive background by just \
        selecting random pairs of genes.
    """
    all_genes = data.var_names.values.astype(str)
    rand_pairs = np.array(list(gen_rand_pairs(all_genes, all_genes, n_pairs)))

    neighbours = List()
    for i in range(data.obsm['spot_neighbours'].shape[0]):
        neighs = np.array(
            data.obsm['spot_neighbours'].values[i, :][0].split(','))
        neighs = neighs[neighs != ''].astype(int)
        neighbours.append(neighs)

    het_vals = np.array([1] * data.shape[0], np.int)

    background = get_lrs_scores(data, rand_pairs, neighbours, het_vals, 0,
                                filter_pairs=False)

    if get_stats:
        pvals, padjs = get_stats(data, background)

        return pvals, padjs, background, rand_pairs

    else:
        return background, rand_pairs

def get_stats(data, background):
    """Gets the stastics based on the background."""
    ###### Generating the stats for the constant background ######
    summ_df = data.uns['lr_summary']
    lr_spot_indices = {}
    spots_to_lr_indices = []
    for i in range(data.shape[0]):
        spots_to_lr_indices.append([])
        for j, lr_ in enumerate(summ_df.index.values):
            spot_indices = np.where(data.obsm['lr_scores'][:,j]>0)[0]
            lr_spot_indices[lr_] = spot_indices
            if i in spot_indices:
                spots_to_lr_indices[i].append(j)

    pvals = np.ones((data.shape[0], summ_df.shape[0]))
    for i, lr_ in enumerate(summ_df.index.values):
        spot_indices = lr_spot_indices[lr_]
        lri_scores = data.obsm['lr_scores'][:, i]
        for j in spot_indices:
            n_greater = len(np.where(background[j, :] >= lri_scores[j])[0])
            n_greater = n_greater if n_greater != 0 else 1  # pseudocount
            pvals[j, i] = n_greater / background.shape[1]

    padjs = np.ones(pvals.shape)
    for spot_i in range(data.shape[0]):
        lr_indices = spots_to_lr_indices[spot_i]
        if len(lr_indices) > 1:
            padjs[spot_i, lr_indices] = multipletests(pvals[spot_i, lr_indices],
                                                      method='fdr_bh')[1]

    return pvals, padjs

def get_rates(lr_scores_, padjs_st_, padjs_const_, total_spots_per_pair, cutoff):
    """Getting proportion of values below different cutoffs"""
    rates_st = []
    rates_const = []
    for j in range(lr_scores_.shape[1]):
        padjs_st_j = padjs_st_[:, j]
        total = total_spots_per_pair[j]
        rates_st.append(len(np.where(padjs_st_j < cutoff)[0]) / total)

        padjs_const_j = padjs_const_[:, j]
        rates_const.append(len(np.where(padjs_const_j < cutoff)[0]) / total)

    rates_st = np.array(rates_st)
    rates_const = np.array(rates_const)
    min_ = int(round(len(rates_st)*.1))
    rates_st = rates_st[np.argsort(-rates_st)[min_:]]
    rates_const = rates_const[np.argsort(-rates_const)[min_:]]
    return np.mean(rates_st), np.mean(rates_const)

def get_fpr_tpr_byPair(data_neg, const_bg_stats, data_pos, const_bg_stats_pos):
    """ Gets the true positive rate & false positive rate by taking the mean of
        the proportion of spots called as significant for a pair in the case of
        negative-control data & true LR pairs.
    """
    # Negative data #
    lr_scores_ = data_neg.obsm['lr_scores']
    padjs_st_ = data_neg.obsm['p_vals']
    padjs_const_ = const_bg_stats['pvals']

    # LR data #
    lr_scores_pos_ = data_pos.obsm['lr_scores']
    padjs_st_pos_ = data_pos.obsm['p_vals']
    padjs_const_pos_ = const_bg_stats_pos['pvals']

    ##### Pval_cutoffs to use for creating ROC curve ######
    n_ = 1000
    pval_cutoffs = np.array(list(range(0, n_))) / n_

    ##### Calculating FPR #######
    total_spots_per_pair = (lr_scores_ > 0).sum(axis=0)
    total_spots_pos_per_pair = (lr_scores_pos_ > 0).sum(axis=0)
    fpr_st_per_pair = []
    fpr_const_per_pair = []
    tpr_st_per_pair = []
    tpr_const_per_pair = []
    for i, cutoff in enumerate(pval_cutoffs):
        # Rates for each of the negative pairs #
        rates_st, rates_const = get_rates(lr_scores_, padjs_st_,
                                              padjs_const_,
                                              total_spots_per_pair, cutoff)
        fpr_st_per_pair.append(rates_st)
        fpr_const_per_pair.append(rates_const)

        # Rates for each of the real pairs #
        rates_st, rates_const = get_rates(lr_scores_pos_, padjs_st_pos_,
                                              padjs_const_pos_,
                                              total_spots_pos_per_pair, cutoff)
        tpr_st_per_pair.append(rates_st)
        tpr_const_per_pair.append(rates_const)

    return fpr_st_per_pair, fpr_const_per_pair, tpr_st_per_pair, tpr_const_per_pair

def diagnostic_scatters(y, group_df, height, fp, lrs, lr_colors,
                        adj, out_plots, file_name):
    """ Diagnostics of relationship between LR expression level & number of
        significant spots.
    """
    x_names = ['median_rank', 'prop_rank']
    colors = vhs.getColors(x_names)
    fig, axes = plt.subplots(ncols=len(x_names),
                             figsize=(4.2 * len(x_names), height))
    axes = axes.ravel()
    for i, x_name in enumerate(x_names):
        ax = axes[i]
        x = group_df.loc[:, x_name].values
        # corr_s = spearmanr(x, y)[0]
        # corr_p = pearsonr(x, y)[0]

        ax.plot(x, y, 'o', c=colors[x_name], )
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlabel(f'LR_{x_name}', fp)
        ax.set_ylabel(f'LR_n_spots_sig', fp)
        # Adding in the LRs #
        for i, lr in enumerate(lrs):
            lr_index = np.where(group_df.index.values == lr)[0][0]
            ax.plot(x[lr_index], y[lr_index], 'o', c=lr_colors[i], )
            ax.text(x[lr_index] + adj[i], y[lr_index],
                    lr.replace('_','-'), c=lr_colors[i])

    fig.suptitle('Relationship between LR expression level & n_sig_spots')
    vhs.dealWithPlot(True, True, True,
                     out_plots, file_name, #f'stlearnBg_LRExprBias_scatter.pdf',
                     300)



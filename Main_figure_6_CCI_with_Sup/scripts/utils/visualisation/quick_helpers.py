""" Helper functions for the quick_plots.py
"""


import numpy as np
import pandas as pd

def get_upset_df(obj_lists, group_names):
    """ Creates the necessary input to draw an upset plot for visualising overlaps
        between multiple groups!
    Args:
        obj_lists (list<list<object>>): List of items in different groups want \
                                          to generate upset-plot for to compare.
        group_names (list<str>): List of strings indicating the names for the \
                                                               different groups.
    Returns:
        pd.DataFrame: This is a dataframe formatted in the required formatted \
                    for input to upsetplot so can visualise multi-overlaps.
    """
    all_hs_genes = []
    [all_hs_genes.extend(de_hs_) for de_hs_ in obj_lists]
    all_hs_genes = np.unique(all_hs_genes)

    de_hs_genes = obj_lists
    samples = group_names

    de_hs_vals = np.zeros((len(samples), len(all_hs_genes)))
    for i, samp in enumerate(samples):
        for j, gene in enumerate(all_hs_genes):
            if gene in de_hs_genes[i]:
                de_hs_vals[i, j] = 1
    de_hs_df = pd.DataFrame(de_hs_vals.transpose(),
                            index=all_hs_genes, columns=samples)

    upset_df = pd.DataFrame()
    col_names = samples
    for idx, col in enumerate(de_hs_df[samples]):
        temp = []
        for i in de_hs_df[col]:
            if i != 0:
                temp.append(True)
            else:
                temp.append(False)
        upset_df[col_names[idx]] = temp

    upset_df['c'] = 1
    example = upset_df.groupby(col_names).count().sort_values('c')

    return example









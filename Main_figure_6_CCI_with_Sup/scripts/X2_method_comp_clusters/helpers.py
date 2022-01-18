"""
Helper functions for comparing the cci results from the clustering.
"""

import numpy as np
import matplotlib.pyplot as plt
from stlearn.plotting.utils import get_colors
import stlearn.plotting.cci_plot_helpers as cph

def chord_plot(int_df, adata, use_label, min_ints: int=2, n_top_ccis: int=10,
				  cmap: str='default', show: bool=True,
				  sig_interactions: bool=True, title=None, label_size: int=10,
				  ):
    """ Chord plot based on interaction df.
    """
    int_df = int_df.transpose()
    fig = plt.figure(figsize=(8, 8))

    flux = int_df.values
    total_ints = flux.sum(axis=1) + flux.sum(axis=0) - flux.diagonal()
    keep = total_ints > min_ints
    # Limit of 10 for good display #
    if sum(keep) > n_top_ccis:
        keep = np.argsort(-total_ints)[0:n_top_ccis]
    flux = flux[:, keep]
    flux = flux[keep, :].astype(float)
    cell_names = int_df.index.values.astype(str)[keep]
    #print(flux[cell_names=='1', cell_names=='5'])
    # Add pseudocount to row/column which has all zeros for the incoming
    # so can make the connection between the two
    # for i in range(flux.shape[0]):
    # 	if np.all(flux[i,:]==0):
    # 		flux[i,flux[:,i]>0] += sys.float_info.min
    # 	elif np.all(flux[:,i]==0):
    # 		flux[flux[i, :] > 0, i] += sys.float_info.min

    #print(cell_names)
    nodes = cell_names

    # Retrieving colors of cell types #
    colors = get_colors(adata, use_label, cmap=cmap, label_set=cell_names)

    ax = plt.axes([0, 0, 1, 1])
    nodePos = cph.chordDiagram(flux, ax, lim=1.25, colors=colors)
    ax.axis('off')
    prop = dict(fontsize=label_size, ha='center', va='center')
    for i in range(len(cell_names)):
        x, y = nodePos[i][0:2]
        ax.text(x, y, nodes[i],
                rotation=nodePos[i][2], size=10, **prop)
    fig.suptitle(title, fontsize=12, fontweight='bold')
    if show:
        plt.show()
    else:
        return fig, ax






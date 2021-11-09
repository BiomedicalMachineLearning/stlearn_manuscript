""" For short, simple plots used for generating diagnostics. """

import numbers
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import upsetplot
import matplotlib
import matplotlib.pyplot as plt

import beautifulcells.visualisation.helpers as vhs

import beautifulcells.visualisation.quick_helpers as qhs

fp = {'weight': 'bold', 'size': 12}

def setUp(font_size=10, axis_font_size=12, fig_title_size=15):
    """Creates bold headings & other matplotlib formatting to create nice plots.
    """
    matplotlib.rcParams.update({'font.size':font_size, 'font.weight': 'bold',
                                'figure.titlesize': fig_title_size,
                                'figure.titleweight': 'bold'})
    global fp
    fp = {'weight': 'bold', 'size': axis_font_size}

def density_scatter(x, y, fig_title='', x_label='', y_label='', figsize=(6.4,4.8),
                    file_name=None, return_ax=False):
    """ Scatter plot with relationship of X & Y, with points coloured by density.
    """
    stack = np.vstack([x, y])
    densities = gaussian_kde(stack)(stack)

    fig, ax = plt.subplots(figsize=figsize)
    fig.suptitle(fig_title)
    ax.scatter(x, y, c=densities)
    ax.set_xlabel(x_label, fp)
    ax.set_ylabel(y_label, fp)
    if not return_ax:
        vhs.dealWithPlot(type(file_name)!=type(None), True, True,
                         '',file_name, 300)
    else:
        return ax

def distrib(x, bins=100, x_label='', fig_title='', log=False, density=False,
            figsize=(6.4,4.8), file_name=None, add_mean=False, logbase=np.e,
            color='blue', alpha=1, ax=None, fig=None, show=True,
            label='', total=None, return_total=False):
    """Plots a histogram of values."""
    if type(ax)==type(None) or type(fig)==type(None):
        fig, ax = plt.subplots(figsize=figsize)

    # Getting the counts in desired format #
    counts, bins = np.histogram(x, bins=bins)
    logcounts = np.log(counts+1)/np.log(logbase) if log else counts
    if density and type(total)==type(None):
        total = sum(logcounts)
        logcounts = logcounts/total
    elif density:
        logcounts = logcounts/total

    ax.hist(bins[:-1], bins, weights=logcounts, color=color, alpha=alpha,
            label=label)
    ax.set_xlabel(x_label, fp)
    if not density:
        ax.set_ylabel(f'log{round(logbase, 2)}-counts' if log else 'counts', fp)
    else:
        ax.set_ylabel('density-'+f'log{round(logbase, 2)}-counts'
                                                       if log else 'counts', fp)
    fig.suptitle(fig_title)

    if add_mean:
        mean = np.mean(x)
        y = ax.get_ylim()[1]*.5
        ax.vlines(mean, 0, y, colors='r')
        ax.text(mean, y, f'mean:{round(mean, 4)}', c='red')

    if show:
        vhs.dealWithPlot(type(file_name)!=type(None), True, True, '',
                                                                 file_name, 300)
    elif not return_total:
        return fig, ax
    else:
        return fig, ax, total

def lineplot(xs, ys, y_names=None, y_colors=None,
             fig_title='', x_label='', y_label='', figsize=(6.4,4.8),
                       file_name=None, return_ax=False, legend_loc='best'):
    """Plots line of multiple data points."""

    if type(y_names)==type(None):
        y_names = list(range(len(ys)))
    if type(y_colors)==type(None):
        y_colors = vhs.getColors(y_names)

    if type(xs[0])!=list and type(xs[0])!=np.array:
        xs = [xs]*len(ys)

    fig, ax = plt.subplots(figsize=figsize)
    for i, y in enumerate(ys):
        ax.plot(xs[i], y, '-', c=y_colors[y_names[i]], linewidth=4,
                label=y_names[i] if type(y_names[i])!=int else None)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xlabel(x_label, fp)
    ax.set_ylabel(y_label, fp)
    fig.suptitle(fig_title)
    plt.legend(loc=legend_loc)

    if not return_ax:
        vhs.dealWithPlot(type(file_name) != type(None), True, True, '',
                         file_name, 300)
    else:
        return ax

def multi_density(values, labels, label_colors=None, bins=100, alpha=.5,
                  fig_title='', x_label='', y_label='', figsize=(6.4,4.8),
                  file_name=None, return_ax=False, legend_loc='best',
                  max_val=None):
    """Plots density of multiple input value across range for comparison."""
    if type(label_colors)==type(None):
        label_colors = vhs.getColors(labels)
    if type(max_val) != type(None):
        values = [np.array(vals)[np.array(vals<max_val)] for vals in values]

    fig, ax = plt.subplots(figsize=figsize)
    for i in range(len(values)):
        ax.hist(values[i], bins, alpha=alpha, label=labels[i], density=True,
                color=label_colors[labels[i]])

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xlabel(x_label, fp)
    ax.set_ylabel(y_label, fp)
    fig.suptitle(fig_title)
    plt.legend(loc=legend_loc)

    if not return_ax:
        vhs.dealWithPlot(type(file_name) != type(None), True, True, '',
                         file_name, 300)
    else:
        return ax

def upset_plot(obj_lists, group_names=None, fig_title='', file_name=None,
               sort_by="cardinality", sort_groups_by=None):
    """ Creates an upset plot, a visualisation which is useful for comparing \
        overlaps between multiple groups when have more than one group.

    Args:
        obj_lists (list<list<object>>): List of items in different groups want \
                                          to generate upset-plot for to compare.
        group_names (list<str>): List of strings indicating the names for the \
                                                               different groups.
    """
    obj_lists = obj_lists[::-1]
    if type(group_names)==type(None):
        group_names = [f'group_{i}' for i in range(len(obj_lists))]
    else:
        group_names = group_names[::-1]

    upset_df = qhs.get_upset_df(obj_lists, group_names)

    upsetplot.plot(upset_df['c'], sort_by=sort_by,
                   sort_categories_by=sort_groups_by)
    plt.title(fig_title, loc='left')
    vhs.dealWithPlot(type(file_name) != type(None), True, True, '',
                                                                 file_name, 300)









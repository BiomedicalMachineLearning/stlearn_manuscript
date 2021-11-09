""" General helper functions common across visualisations.
"""

import matplotlib
import matplotlib.pyplot as plt

import numpy
import collections

import beautifulcells.utils.inputHelpers as input
import beautifulcells.utils.calculations as calcs

def dealWithPlot(savePlot, showPlot, closePlot, folder, plotName, dpi,
				 tightLayout=True):
	""" Deals with the current matplotlib.pyplot.
	"""

	if tightLayout:
		plt.tight_layout()

	if savePlot:
		plt.savefig(folder+plotName, dpi=dpi,
					format=plotName.split('.')[-1])

	if showPlot:
		plt.show()

	if closePlot:
		plt.close()

def getFigAxes(fig, axes, figsize, **kwargs):
	if type(fig) == type(None):
		fig = plt.figure(figsize=figsize, **kwargs)
	if type(axes) == type(None):
		axes = fig.add_subplot(111)  # 1 row, 1 col, start at cell 1.

	return fig, axes

def getColors(labels, labelSet=None, colorMap='tab20', rgb=False):
	""" Gets an OrderedDict of colors; the order indicates the frequency of \
	labels from largest to smallest.

	Args:
		labels (numpy.array<str>): Indicates a set of labels for observations.

		labelSet (list-like<str>): Indicates the set of labels in labels. \
									If None, calculated based on labels.

		colorMap (str): A matplotlib colormap.

		rgb (bool): If True, colors indicated by rgb value, if false hexcode.

	Returns:
		dict<str, tuple or str>: An ordered dict indicating the labels which \
					occur most to least frequently and the associated colors.
	"""
	# Determining the set of labels #
	labelSet = input.returnDefaultIfNone(labelSet,
										 calcs.getOrderedLabelSet(labels))

	# Initialising the ordered dict #
	cellTypeColors = {}

	# Ordering the cells according to their frequency and obtaining colors #
	nLabels = len(labelSet)
	cmap = plt.cm.get_cmap(colorMap, nLabels)
	rgbs = [cmap(i)[:3] for i in range(nLabels)]
	#rgbs = list(numpy.array(rgbs)[order]) # Make sure color order is the same.

	# Populating the color dictionary with rgb values or hexcodes #
	for i in range(len(labelSet)):
		cellType = labelSet[i]
		rgbi = rgbs[i]
		if not rgb:
			cellTypeColors[cellType] = matplotlib.colors.rgb2hex(rgbi)
		else:
			cellTypeColors[cellType] = rgbi

	return cellTypeColors


#!/usr/bin/env python
"""
Plot mRNA counts

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/27/2014
"""

from __future__ import division

import argparse
import os

import tables
import numpy as np
from scipy import stats
import matplotlib
# matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.utils.knowledgebase_fixture_manager import loadKnowledgeBase

# TODO: account for complexation

def main(simOutDir, plotOutDir, plotOutFileName):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get the names of rnas from the KB

	kb = loadKnowledgeBase(simOutDir)

	rnaIds = kb.rnaData["id"][kb.rnaIndexToMonomerMapping]

	proteinIds = kb.monomerData["id"]

	with tables.open_file(os.path.join(simOutDir, "BulkMolecules.hdf")) as bulkMoleculesFile:

		names = bulkMoleculesFile.root.names
		bulkMolecules = bulkMoleculesFile.root.BulkMolecules

		moleculeIds = names.moleculeIDs.read()

		rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)

		rnaCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, rnaIndexes]

		proteinIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in proteinIds], np.int)

		proteinCountsBulk = bulkMolecules.read(0, None, 1, "counts")[:, proteinIndexes]

	relativeMRnaCounts = rnaCountsBulk[-1, :] #/ rnaCountsBulk[-1, :].sum()
	relativeProteinCounts = proteinCountsBulk[-1, :] #/ proteinCountsBulk[-1, :].sum()

	plt.plot(relativeMRnaCounts, relativeProteinCounts, '.')

	# scatterPlotWithHistograms(relativeMRnaCounts, relativeProteinCounts)

	plt.show()

	# plt.plot(rnaCountsBulk.flatten(), proteinCountsBulk.flatten(), '.')

	# plt.show()

	import ipdb; ipdb.set_trace()

	# plt.figure(figsize = (8.5, 11))

	# plt.plot(time / 60, nActive)
	# plt.plot([time[0] / 60., time[-1] / 60.], [2 * nActive[0], 2 * nActive[0]], "r--")
	# plt.xlabel("Time (min)")
	# plt.ylabel("Counts")
	# plt.title("Active Ribosomes Final:Initial = %0.2f" % (nActive[-1] / float(nActive[0])))

	# plt.savefig(os.path.join(plotOutDir, plotOutFileName))


def scatterPlotWithHistograms(x, y):
	# Adapted directly from
	# http://matplotlib.org/examples/axes_grid/scatter_hist.html

	# import numpy as np
	# import matplotlib.pyplot as plt
	from mpl_toolkits.axes_grid1 import make_axes_locatable

	# # the random data
	# x = np.random.randn(1000)
	# y = np.random.randn(1000)


	fig, axScatter = plt.subplots(figsize=(5.5,5.5))

	# the scatter plot:
	axScatter.scatter(x, y)
	axScatter.set_aspect(1.)

	# create new axes on the right and on the top of the current axes
	# The first argument of the new_vertical(new_horizontal) method is
	# the height (width) of the axes to be created in inches.
	divider = make_axes_locatable(axScatter)
	axHistx = divider.append_axes("top", 1.2, pad=0.1, sharex=axScatter)
	axHisty = divider.append_axes("right", 1.2, pad=0.1, sharey=axScatter)

	# make some labels invisible
	plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
	         visible=False)

	# now determine nice limits by hand:
	binwidth = 0.25
	xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
	lim = ( int(xymax/binwidth) + 1) * binwidth

	bins = np.arange(-lim, lim + binwidth, binwidth)
	axHistx.hist(x, bins=bins)
	axHisty.hist(y, bins=bins, orientation='horizontal')

	# the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
	# thus there is no need to manually adjust the xlim and ylim of these
	# axis.

	#axHistx.axis["bottom"].major_ticklabels.set_visible(False)
	for tl in axHistx.get_xticklabels():
	    tl.set_visible(False)
	axHistx.set_yticks([0, 50, 100])

	#axHisty.axis["left"].major_ticklabels.set_visible(False)
	for tl in axHisty.get_yticklabels():
	    tl.set_visible(False)
	axHisty.set_xticks([0, 50, 100])

	plt.draw()
	# plt.show()


if __name__ == "__main__":
	# parser = argparse.ArgumentParser()
	# parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	# parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	# parser.add_argument("plotOutFileName", help = "File name to produce", type = str)

	# args = parser.parse_args().__dict__

	# main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"])

	# TODO: reinstate this analysis
	pass

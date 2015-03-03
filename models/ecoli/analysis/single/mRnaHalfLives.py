#!/usr/bin/env python
"""
Plot mRNA half lives (observed vs. actual)

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/30/2015
"""

from __future__ import division

import argparse
import os

import tables
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

import wholecell.utils.constants
from wholecell.io.tablereader import TableReader

# TODO: account for complexation

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get the names of rnas from the KB

	kb = cPickle.load(open(kbFile, "rb"))

	isMRna = kb.rnaData["isMRna"]
	rnaIds = kb.rnaData["id"][isMRna]

	rnaDegradationListenerFile = TableReader(os.path.join(simOutDir, "RnaDegradationListener"))
	time = rnaDegradationListenerFile.readColumn("time")
	countRnaDegraded = rnaDegradationListenerFile.readColumn('countRnaDegraded')
	rnaDegradationListenerFile.close()

	plt.figure(figsize = (8.5, 11))

	# Assuming that decay rates are the average of instantaneous degradation rates (simulation time points)
	rnaDegradationRate = countRnaDegraded[1:,:].sum(axis = 0)[isMRna] / 3600.
	# Assuming that decay rates are the median of instantaneous degradation rates (simulation time points)
	# rnaDegradationRate = countRnaDegraded[1:,:].median(axis = 0)[isMRna]

	rnaDegradationRateStd = countRnaDegraded[1:,:].std(axis = 0)[isMRna]


	expectedDegradationRate = kb.rnaData['degRate'][isMRna].asNumber()

	maxLine = 1.1 * max(expectedDegradationRate.max(), rnaDegradationRate.max())
	
	plt.plot([0, maxLine], [0, maxLine], '--r')	
	plt.plot(expectedDegradationRate, rnaDegradationRate, 'o', markeredgecolor = 'k', markerfacecolor = 'none')
	#plt.errorbar(expectedDegradationRate, rnaDegradationRate, rnaDegradationRateStd)
	Correlation_ExpPred = np.corrcoef(expectedDegradationRate, rnaDegradationRate)[0][1]
	print "Correlation expected and predicted half-lives = %.3f" % Correlation_ExpPred

	plt.xlabel("Expected RNA decay")
	plt.ylabel("Actual RNA decay (at final time step)")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)


if __name__ == "__main__":
	defaultKBFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--kbFile", help = "KB file name", type = str, default = defaultKBFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["kbFile"])

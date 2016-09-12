#!/usr/bin/env python
"""
Plot TF probabilities

@author: Travis Horst
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/22/16
"""

from __future__ import division

import argparse
import os

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import cPickle

import wholecell.utils.constants
from wholecell.io.tablereader import TableReader

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	sim_data = cPickle.load(open(simDataFile, "rb"))

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	
	tfProb = TableReader(os.path.join(simOutDir, "RnaSynthProb"))
	pTfBound = tfProb.readColumn("pTfBound")
	pPromoterBound = tfProb.readColumn("pPromoterBound")
	nTfBound = tfProb.readColumn("nTfBound")
	nPromoterBound = tfProb.readColumn("nPromoterBound")
	nActualBound = tfProb.readColumn("nActualBound")
	tfProb.close()

	trpActiveExpectedProbBound = sim_data.pPromoterBound["CPLX-125__active"]["CPLX-125"]
	trpInactiveExpectedProbBound = sim_data.pPromoterBound["CPLX-125__inactive"]["CPLX-125"]
	tyrActiveExpectedProbBound = sim_data.pPromoterBound["MONOMER0-162__active"]["MONOMER0-162"]
	tyrInactiveExpectedProbBound = sim_data.pPromoterBound["MONOMER0-162__inactive"]["MONOMER0-162"]

	recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames
	tfs = sorted(set([x.split("__")[-1] for x in recruitmentColNames if x.split("__")[-1] != "alpha"]))
	trpRIndex = [i for i, tf in enumerate(tfs) if tf == "CPLX-125"][0]
	tyrRIndex = [i for i, tf in enumerate(tfs) if tf == "MONOMER0-162"][0]

	plt.figure(figsize = (8.5, 11))
	ax = plt.subplot(2,2,1)
	plt.plot(time, np.array([pTfBound[:,trpRIndex], pPromoterBound[:,trpRIndex]]).T)
	plt.plot([0, time[-1]], [trpActiveExpectedProbBound, trpActiveExpectedProbBound])
	plt.plot([0, time[-1]], [trpInactiveExpectedProbBound, trpInactiveExpectedProbBound])
	plt.ylabel("TrpR Probability Bound")
	plt.legend(["TF", "Promoter", "Expected Promoter - Active", "Expected Promoter - Inactive"], fontsize = 6)
	ax.set_ylim([0, 1.2])

	ax = plt.subplot(2,2,2)
	plt.plot(time, np.array([pTfBound[:,tyrRIndex], pPromoterBound[:,tyrRIndex]]).T)
	plt.plot([0, time[-1]], [tyrActiveExpectedProbBound, tyrActiveExpectedProbBound])
	plt.plot([0, time[-1]], [tyrInactiveExpectedProbBound, tyrInactiveExpectedProbBound])
	plt.ylabel("TyrR Probability Bound")
	plt.legend(["TF", "Promoter", "Expected Promoter - Active", "Expected Promoter - Inactive"], fontsize = 6)
	ax.set_ylim([0, 1.2])

	ax = plt.subplot(2,2,3)
	plt.plot(time, np.array([nActualBound[:,trpRIndex], nTfBound[:,trpRIndex]]).T)
	plt.ylabel("TrpR Number Bound")
	plt.legend(["Actual", "Expected"], fontsize = 10)
	ax.set_ylim([0, max(nActualBound[:,trpRIndex]) + 2])

	ax = plt.subplot(2,2,4)
	plt.plot(time, np.array([nActualBound[:,tyrRIndex], nTfBound[:,tyrRIndex]]).T)
	plt.ylabel("TyrR Number Bound")
	plt.legend(["Actual", "Expected"], fontsize = 10)
	ax.set_ylim([0, max(nActualBound[:,tyrRIndex]) + 2])

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)

if __name__ == "__main__":
	defaultSimDataFile = os.path.join(
			wholecell.utils.constants.SERIALIZED_KB_DIR,
			wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME
			)

	parser = argparse.ArgumentParser()
	parser.add_argument("simOutDir", help = "Directory containing simulation output", type = str)
	parser.add_argument("plotOutDir", help = "Directory containing plot output (will get created if necessary)", type = str)
	parser.add_argument("plotOutFileName", help = "File name to produce", type = str)
	parser.add_argument("--simDataFile", help = "KB file name", type = str, default = defaultSimDataFile)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])

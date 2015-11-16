#!/usr/bin/env python
"""
Plot reaction max rate over course of the simulation.

@date: Created 7/02/2015
@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import division

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, kbFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	enzymeKineticsdata = TableReader(os.path.join(simOutDir, "EnzymeKinetics"))
	
	enzymeKineticsArray = enzymeKineticsdata.readColumn("reactionRates")
	perEnzymeRates = enzymeKineticsdata.readColumn("perEnzymeRates")

	reactionIDs = enzymeKineticsdata.readAttribute("reactionIDs")
	
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	
	enzymeKineticsdata.close()

	reactionRateArray = np.transpose(enzymeKineticsArray)
	perEnzymeRateArray = np.transpose(perEnzymeRates)

	plt.figure(figsize = (8.5, 11))
	# Reaction rate
	plt.subplot(2,1,1)

	plt.title("Enzyme Kinetics")

	lineLabels = []

	i = 0
	for timeCourse in reactionRateArray:
		if (np.amax(timeCourse) < np.inf) and (i < len(reactionIDs)):
			plt.plot(time / 60, timeCourse)
			lineLabels.append(reactionIDs[i][:15])
		i += 1

	plt.xlabel("Time (min)")
	plt.ylabel("Reaction Rate (reactions/second)")
	plt.legend(lineLabels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	# Per-enzyme reaction rate
	
	lineLabels = []
	i = 0
	plt.subplot(2,1,2)
	for timeCourse in perEnzymeRateArray:
		if (np.amax(timeCourse) < np.inf) and (i < len(reactionIDs)):
			plt.plot(time / 60, timeCourse)
			lineLabels.append(reactionIDs[i][:15])
		i += 1

	plt.xlabel("Time (min)")
	plt.ylabel("Per Enzyme Rate (reactions/enzyme-second)")

	plt.subplots_adjust(left=.12, right=.65, top=0.9, bottom=0.1)

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")

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

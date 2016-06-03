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

from models.ecoli.processes.metabolism import COUNTS_UNITS, VOLUME_UNITS, TIME_UNITS

COLORS = [
	[166,206,227],
	[31,120,180],
	[178,223,138],
	[51,160,44],
	[251,154,153],
	[227,26,28],
	[253,191,111],
	[255,127,0],
	[202,178,214],
	[106,61,154],
	[255,255,153],
	[177,89,40]]

CMAP_COLORS = [[shade/255. for shade in color] for color in COLORS]

MAX_STRLEN = 30

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	dualValues = np.array(fbaResults.readColumn("dualValues")).T
	moleculeIDs = np.array(fbaResults.readAttribute("outputMoleculeIDs"))
	fbaResults.close()

	num_values = 3

	plt.figure(figsize = (8.5, 11))
	plt.title("FBA Dual Values")

	# Get the num_values highest reduced price values at each timestep
	highest_partition = np.argpartition(dualValues,-num_values,axis=0)[-num_values:].T
	for plotNum, molIdx in enumerate(np.unique(highest_partition)):
		plt.plot(time / 60., dualValues[molIdx], '--', color=CMAP_COLORS[plotNum%len(CMAP_COLORS)], label=moleculeIDs[molIdx][:MAX_STRLEN])

	# num_values lowest reduced price values at each timestep
	lowest_partition = np.argpartition(dualValues,num_values,axis=0)[:num_values].T
	for plotNum, molIdx in enumerate(np.unique(lowest_partition)):
		plt.plot(time / 60., dualValues[molIdx], '.', color=CMAP_COLORS[plotNum%len(CMAP_COLORS)], label=moleculeIDs[molIdx][:MAX_STRLEN])

	plt.xlabel("Time (min)")
	plt.ylabel("Reduced Cost (All molecules appearing in the top or bottom {} reduced cost for a timestep)".format(num_values))
	plt.legend(framealpha=.5, fontsize=6, loc='best')

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName)
	plt.close("all")

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
#!/usr/bin/env python
"""
Plot tRNA counts

@author: Heejo Choi
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/19/2017
"""

import argparse
import os
import cPickle
import numpy as np
import matplotlib.pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):
	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Get time
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")

	# Get tRNA IDs and counts
	sim_data = cPickle.load(open(simDataFile, "rb"))
	isTRna = sim_data.process.transcription.rnaData["isTRna"]
	rnaIds = sim_data.process.transcription.rnaData["id"][isTRna]
	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))
	moleculeIds = bulkMolecules.readAttribute("objectNames")
	rnaIndexes = np.array([moleculeIds.index(moleculeId) for moleculeId in rnaIds], np.int)
	rnaCountsBulk = bulkMolecules.readColumn("counts")[:, rnaIndexes]
	bulkMolecules.close()

	# Plot
	fig = plt.figure(figsize = (8.5, 11))
	ax = plt.subplot(1, 1, 1)
	ax.plot(time, rnaCountsBulk)
	ax.set_xlim([time[0], time[-1]])
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Counts of tRNAs")
	ax.spines["right"].set_visible(False)
	ax.spines["top"].set_visible(False)
	ax.tick_params(right = "off", top = "off", which = "both", direction = "out")

	from wholecell.analysis.analysis_tools import exportFigure
	exportFigure(plt, plotOutDir, plotOutFileName, metadata)
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
	parser.add_argument("--validationDataFile", help = "Validation data file name", type = str)

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])

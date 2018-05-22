#!/usr/bin/env python
"""
Plot enzymatic capacity of tryptophan synthase vs amount of tryptophan needed by translation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/22/2017
"""

import argparse
import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants

BURN_IN = 10

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	sim_data = cPickle.load(open(simDataFile, "r"))
	trpIdx = sim_data.moleculeGroups.aaIDs.index("TRP[c]")

	growthLimits = TableReader(os.path.join(simOutDir, "GrowthLimits"))

	trpRequests = growthLimits.readColumn("aaRequestSize")[BURN_IN:, trpIdx]

	growthLimits.close()

	bulkMolecules = TableReader(os.path.join(simOutDir, "BulkMolecules"))

	moleculeIds = bulkMolecules.readAttribute("objectNames")

	trpSynIdx = moleculeIds.index("TRYPSYN[c]")

	trpSynCounts = bulkMolecules.readColumn("counts")[BURN_IN:, trpSynIdx]

	bulkMolecules.close()

	trpSynKcat = 2**( (37. - 25.) / 10.) * 4.1 # From PMID 6402362 (kcat of 4.1/s measured at 25 C)

	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time")[BURN_IN:] - initialTime
	timeStep = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")[BURN_IN:]


	trpSynMaxCapacity = trpSynKcat * trpSynCounts * timeStep

	plt.figure(figsize = (8.5, 11))

	plt.subplot(3, 1, 1)

	plt.plot(time / 60., trpSynMaxCapacity, linewidth = 2)
	plt.ylabel("Tryptophan Synthase Max Capacity")

	plt.subplot(3, 1, 2)

	plt.plot(time / 60., trpRequests, linewidth = 2)
	plt.ylabel("TRP requested by translation")

	plt.subplot(3, 1, 3)

	plt.plot(time / 60., trpSynMaxCapacity / trpRequests, linewidth = 2)
	plt.xlabel("Time (min)")
	plt.ylabel("(Max capacity) / (Request)")


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

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"])

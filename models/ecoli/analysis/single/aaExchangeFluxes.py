#!/usr/bin/env python

import argparse
import os
import cPickle

import numpy as np
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata=None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	# Amino acid IDs
	sim_data = cPickle.load(open(simDataFile, "rb"))
	aaIDs = sim_data.moleculeGroups.aaIDs

	# Amino acid exchanges fluxes
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
	externalExchangeFluxes = fba_results.readColumn("externalExchangeFluxes")
	externalMoleculeIDs = fba_results.readAttribute("externalMoleculeIDs")

	# Plot
	rows = 6
	cols = 4
	fig = plt.figure(figsize = (8, 11.5))

	for plotIndex, aa in enumerate(aaIDs):
		ax = plt.subplot(rows, cols, plotIndex + 1)

		if not aa.startswith("L-SELENOCYSTEINE"):
			aa = aa[:-3] + "[p]"
		if aa in externalMoleculeIDs:
			aaFlux = externalExchangeFluxes[:, externalMoleculeIDs.index(aa)]
		else:
			aaFlux = np.zeros(len(time))

		ax.plot(time / 60., aaFlux)
		ax.set_xlabel("Time (min)", fontsize = 6)
		ax.set_ylabel("mmol/gDCW/hr", fontsize = 6)
		ax.set_title("%s" % aa, fontsize = 6, y = 1.1)
		ax.tick_params(which = "both", direction = "out", labelsize = 6)

	plt.rc("font", size = 6)
	plt.suptitle("External exchange fluxes of amino acids", fontsize = 10)

	plt.subplots_adjust(hspace = 1, wspace = 1)
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
	parser.add_argument("--validationDataFile")

	args = parser.parse_args().__dict__

	main(args["simOutDir"], args["plotOutDir"], args["plotOutFileName"], args["simDataFile"], args["validationDataFile"])

#!/usr/bin/env python

import argparse
import os
import cPickle

import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from wholecell.io.tablereader import TableReader
import wholecell.utils.constants
from wholecell.utils import units

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata=None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)

	# Plot glucose exchange flux
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime

	fba_results = TableReader(os.path.join(simOutDir, "FBAResults"))
	exFlux = fba_results.readColumn("externalExchangeFluxes")
	exMolec = fba_results.readAttribute("externalMoleculeIDs")
	glcFlux = exFlux[:,exMolec.index("GLC[p]")]
	oxygenFlux = exFlux[:,exMolec.index("OXYGEN-MOLECULE[p]")]

	ax1.plot(time / 60. / 60., -1. * glcFlux, label="Glucose exchange flux coefficient")
	ax1.set_ylabel("External\nglucose\n(mmol/gDCW/hr)", fontsize = 8)
	ax1.set_title("GLC[p]", fontsize = 8)
	ax1.tick_params(labelsize =8, which = "both", direction = "out")

	ax2.plot(time / 60. / 60., -1. * oxygenFlux, label="Oxygen exchange flux coefficient")
	ax2.set_title("OXYGEN-MOLECULE[p]", fontsize = 8)
	ax2.set_ylabel("External\noxygen\n(mmol/gDCW/hr)", fontsize = 8)
	ax2.set_xlabel("Time (hr)", fontsize = 8)
	ax2.tick_params(labelsize = 8, which = "both", direction = "out")


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

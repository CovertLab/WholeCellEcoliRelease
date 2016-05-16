#!/usr/bin/env python
"""
@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/8/2014
"""

from __future__ import division

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

GLUCOSE_ID = "GLC[p]"

FLUX_UNITS = units.mmol / units.g / units.h
MASS_UNITS = units.fg
GROWTH_UNITS = units.fg / units.s

def main(simOutDir, plotOutDir, plotOutFileName, simDataFile, validationDataFile, metadata = None):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	sim_data = cPickle.load(open(simDataFile, "rb"))

	fbaResults = TableReader(os.path.join(simOutDir, "FBAResults"))
	initialTime = TableReader(os.path.join(simOutDir, "Main")).readAttribute("initialTime")
	time = TableReader(os.path.join(simOutDir, "Main")).readColumn("time") - initialTime
	timeStepSec = TableReader(os.path.join(simOutDir, "Main")).readColumn("timeStepSec")
	externalExchangeFluxes = fbaResults.readColumn("externalExchangeFluxes")

	externalMoleculeIDs = np.array(fbaResults.readAttribute("externalMoleculeIDs"))

	fbaResults.close()

	glucoseIdx = np.where(externalMoleculeIDs == GLUCOSE_ID)[0][0]
	glucoseFlux = FLUX_UNITS * externalExchangeFluxes[:, glucoseIdx]

	mass = TableReader(os.path.join(simOutDir, "Mass"))
	cellMass = MASS_UNITS * mass.readColumn("cellMass")
	cellDryMass = MASS_UNITS * mass.readColumn("dryMass")
	growth = GROWTH_UNITS * mass.readColumn("growth") / timeStepSec
	mass.close()

	glucoseMW = sim_data.getter.getMass([GLUCOSE_ID])[0]

	glucoseMassFlux = glucoseFlux * glucoseMW * cellDryMass

	glucoseMassYield = growth / -glucoseMassFlux

	fig = plt.figure(figsize = (8.5, 11))
	plt.plot(time, glucoseMassYield)
	plt.xlabel("Time (s)")
	plt.ylabel("g cell / g glucose")

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

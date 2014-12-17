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

import tables
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

import wholecell.utils.constants
from wholecell.utils import units

GLUCOSE_ID = "GLC-D[e]"

FLUX_UNITS = units.mmol / units.L / units.s
MASS_UNITS = units.fg
GROWTH_UNITS = units.fg / units.s

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	kb = cPickle.load(open(kbFile, "rb"))

	with tables.open_file(os.path.join(simOutDir, "FBAResults.hdf")) as h5file:
		time = h5file.root.FBAResults.col("time")
		timeStep = h5file.root.FBAResults.col("timeStep")
		externalExchangeFluxes = h5file.root.FBAResults.col("externalExchangeFluxes")

		names = h5file.root.names
		externalMoleculeIDs = np.array(names.externalMoleculeIDs.read())

	glucoseIdx = np.where(externalMoleculeIDs == GLUCOSE_ID)[0][0]
	glucoseFlux = FLUX_UNITS * externalExchangeFluxes[:, glucoseIdx]

	with tables.open_file(os.path.join(simOutDir, "Mass.hdf")) as h5file:
		table = h5file.root.Mass
		cellMass = MASS_UNITS * table.read(0, None, 1, "cellMass")
		cellDryMass = MASS_UNITS * table.read(0, None, 1, "dryMass")
		growth = GROWTH_UNITS * table.read(0, None, 1, "growth")

	cellDensity = kb.cellDensity
	glucoseMW = kb.getMass([GLUCOSE_ID])[0]

	glucoseMassFlux = glucoseFlux * glucoseMW * cellMass / cellDensity

	glucoseMassYield = growth / glucoseMassFlux

	fig = plt.figure(figsize = (8.5, 11))
	plt.plot(time, glucoseMassYield)
	plt.xlabel("Time (s)")
	plt.ylabel("g cell / g glucose")

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

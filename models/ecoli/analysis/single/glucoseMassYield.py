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

FLUX_UNITS = "M/s"

def main(simOutDir, plotOutDir, plotOutFileName, kbFile):

	if not os.path.isdir(simOutDir):
		raise Exception, "simOutDir does not currently exist as a directory"

	if not os.path.exists(plotOutDir):
		os.mkdir(plotOutDir)

	kb = cPickle.load(open(kbFile, "rb"))

	with tables.open_file(os.path.join(simOutDir, "FBAResults.hdf")) as h5file:
		time = h5file.root.FBAResults.col("time")
		timeStep = h5file.root.FBAResults.col("timeStep")
		outputFluxes = h5file.root.FBAResults.col("outputFluxes")

		names = h5file.root.names
		outputMoleculeIDs = np.array(names.outputMoleculeIDs.read())

	glucoseIdx = np.where(outputMoleculeIDs == 'GLU-L[c]')[0][0]
	glucoseFlux = outputFluxes[:,glucoseIdx] #mol/L/s

	with tables.open_file(os.path.join(simOutDir, "Mass.hdf")) as h5file:
		table = h5file.root.Mass
		cellMass = np.array([x["cellMass"] for x in table.iterrows()]) # fg
		cellDryMass = np.array([x["dryMass"] for x in table.iterrows()]) # fg
		growth = np.array([x["growth"] for x in table.iterrows()]) # fg

	cellDensity = kb.cellDensity.asNumber(units.fg/units.L)
	glucoseMW = np.sum(kb.bulkMolecules['mass'][kb.bulkMolecules['moleculeId'] == 'GLU-L[c]']).asNumber(units.g/units.mol)

	glucoseMassFlux = glucoseFlux * glucoseMW * cellDryMass / cellDensity * 10**15 # fg glucose / s

	massGrowth = growth / cellMass # fg / s

	glucoseMassYield = massGrowth / glucoseMassFlux

	fig = plt.figure(figsize = (8.5, 11))
	plt.plot(time, glucoseMassYield)
	plt.xlabel("Time (s)")
	plt.ylabel("g cell / g glucose")

	plt.savefig(os.path.join(plotOutDir, plotOutFileName))


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

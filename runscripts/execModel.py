#!/usr/bin/env python

"""
execModel.py

Executes a model.

"""

from __future__ import division

from wholecell.sim.simulation import getSimOptsFromEnvVars
from models.ecoli.sim.simulation import EcoliSimulation
from models.ecoli_metabolism.sim.simulation import EcoliMetabolismSimulation
import wholecell.utils.constants
import os
import argparse

def main(modelLevel, kbDirectory, simDirectory):

	simOpts = getSimOptsFromEnvVars(
		["outputDir", "kbLocation", "logToDisk",
		"logToShell", "logToDiskEvery", "overwriteExistingFiles"
		]
		)

	fitKbFileName = (
		wholecell.utils.constants.SERIALIZED_KB_PREFIX +
		("_Fit_%d" % modelLevel) +
		wholecell.utils.constants.SERIALIZED_KB_SUFFIX
		)

	simOpts["kbLocation"] = os.path.join(kbDirectory, fitKbFileName)
	simOpts["outputDir"] = os.path.join(simDirectory, "model_level_%d" % modelLevel, "simOut")
	simOpts["logToDisk"] = True
	simOpts["logToShell"] = False
	simOpts["logToDiskEvery"] = 10
	simOpts["overwriteExistingFiles"] = False

	if not os.path.exists(simOpts["outputDir"]):
		os.makedirs(simOpts["outputDir"])

	if modelLevel == 1:
		sim = EcoliMetabolismSimulation(**simOpts)

	if modelLevel == 2:
		simOpts["logToDisk"] = False
		simOpts["logToShell"] = True
		sim = EcoliSimulation(**simOpts)

	sim.run()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument(
		"modelLevel",
		help = "Stage of model execution",
		type = int)
	parser.add_argument(
		"kbDirectory",
		help = "Directory containing kbs",
		type = str)
	parser.add_argument(
		"simDirectory",
		help = "Directory containing simulations",
		type = str)

	args = parser.parse_args().__dict__

	main(args["modelLevel"], args["kbDirectory"], args["simDirectory"])

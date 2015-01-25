#!/usr/bin/env python

"""
execModel.py

Executes a model.

"""

from __future__ import division

import os
import argparse
import shutil

from wholecell.sim.simulation import getSimOptsFromEnvVars
from models.ecoli.sim.simulation import EcoliSimulation
from models.ecoli_metabolism.sim.simulation import EcoliMetabolismSimulation
import wholecell.utils.constants

def main(modelLevel, kbDirectory, simDirectory):
	# TODO: bitch at Derek about all the wacky assumptions in this script

	simOpts = getSimOptsFromEnvVars()

	fitKbFileName = (
		wholecell.utils.constants.SERIALIZED_KB_PREFIX +
		("_Fit_%d" % modelLevel) +
		wholecell.utils.constants.SERIALIZED_KB_SUFFIX
		)

	simOpts["seed"] = 0
	simOpts["kbLocation"] = os.path.join(kbDirectory, fitKbFileName)
	simOpts["outputDir"] = os.path.join(simDirectory, "model_level_%d" % modelLevel, "simOut")
	simOpts["logToDisk"] = True
	simOpts["logToShell"] = False
	simOpts["logToDiskEvery"] = 10
	simOpts["overwriteExistingFiles"] = False

	# Clear out old directory, if it exists
	shutil.rmtree(simOpts["outputDir"], ignore_errors = True)

	metadataDir = os.path.join(simDirectory, "model_level_%d" % modelLevel, "metadata")

	if not os.path.exists(simOpts["outputDir"]):
		os.makedirs(simOpts["outputDir"])

	if not os.path.exists(metadataDir):
		os.makedirs(metadataDir)

	if modelLevel == 1:
		simOpts["lengthSec"] = 3600
		sim = EcoliMetabolismSimulation(**simOpts)
		EcoliMetabolismSimulation.printAnalysisSingleFiles(fileName = os.path.join(metadataDir, "singleAnalysis.list"))

	if modelLevel == 2:
		simOpts["logToDisk"] = False
		simOpts["logToShell"] = True
		sim = EcoliSimulation(**simOpts)
		EcoliSimulation.printAnalysisSingleFiles(fileName = os.path.join(metadataDir, "singleAnalysis.list"))


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

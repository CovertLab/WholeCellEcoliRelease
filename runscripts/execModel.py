#!/usr/bin/env python

"""
execModel.py

Executes a model.

"""

from __future__ import division

from wholecell.sim.simulation import getSimOptsFromEnvVars
from models.ecoli.sim.simulation import EcoliSimulation
import wholecell.utils.constants
import os
import argparse

def main(modelLevel, kbDirectory, outputDirectory):

	simOpts = getSimOptsFromEnvVars(
		["outputDir", "kbLocation"]
		)

	fitKbFileName = (
		wholecell.utils.constants.SERIALIZED_KB_PREFIX +
		("_Fit_%d" % modelLevel) +
		wholecell.utils.constants.SERIALIZED_KB_SUFFIX
		)

	simOpts["kbLocation"] = os.path.join(kbDirectory, fitKbFileName)
	simOpts["outputDir"] = os.path.join(outputDirectory,"model_level_%d" % modelLevel)

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
		"outputDirectory",
		help = "Directory containing output",
		type = str)

	args = parser.parse_args().__dict__

	main(args["modelLevel"], args["kbDirectory"], args["outputDirectory"])

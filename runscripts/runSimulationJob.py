#!/usr/bin/env python

"""
runSimulationJob.py

Runs a simulation, called from a job script.

Options are passed using environmental variables.

"""

from __future__ import division

from wholecell.sim.simulation import getSimOptsFromEnvVars
from models.ecoli.sim.simulation import EcoliSimulation
import argparse
import os
import json

def main(submissionTime, variant, variantId):

	simOpts = getSimOptsFromEnvVars(
			["outputDir", "logToDisk", "overwriteExistingFiles"]
			)

	# Define simulation options specific to this script
	simOpts["outputDir"] = os.path.join(
		"out",
		"%s" % submissionTime,
		"%s_%06d" % (variant, variantId),
		"%06d" % simOpts["seed"],
		"simOut",
		)
	simOpts["logToDisk"] = True
	simOpts["overwriteExistingFiles"] = False

	sim = EcoliSimulation(**simOpts)

	sim.run()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("submissionTime", help = "Time of job submission", type = str)
	parser.add_argument("variant", help = "Variant", type = str)
	parser.add_argument("variantId", help = "Variant id", type = int)


	args = parser.parse_args().__dict__

	main(args["submissionTime"], args["variant"], args["variantId"])
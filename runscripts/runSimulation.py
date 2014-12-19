#!/usr/bin/env python

'''
runSimulation.py

Runs a simulation.

Run a simulation using default parameters:
~/wcEcoli$ python runSimulation.py

Run a simulation using a provided configuration (JSON) file:
~/wcEcoli$ python runSimulation.py simParameters.json

'''

from __future__ import division

from wholecell.sim.simulation import getSimOptsFromEnvVars
from models.ecoli.sim.simulation import EcoliSimulation
import argparse
import os
import json

def main(submissionTime):
	
	simOpts = getSimOptsFromEnvVars(
		["outputDir", "logToDisk", "overwriteExistingFiles"]
		)

	# Define simulation options specific to this script
	simOpts["outputDir"] = os.path.join(
		"out",
		"%s" % submissionTime,
		"wildtype_000000",
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

	args = parser.parse_args().__dict__

	main(args["submissionTime"])

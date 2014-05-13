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

import wholecell.sim.simulation
import argparse
import os

def main(submissionTime):

	seed = 0

	if os.environ.has_key("PBS_ARRAYID"):
		seed = int(os.environ["PBS_ARRAYID"])

	outputDir = os.path.join("out", "simOut", "%s" % submissionTime, "%06d" % seed)

	simOpts = dict(
		seed = seed,
		lengthSec = 10,
		logToShell = True,
		logToDisk = True,
		logToDiskEvery = 10,
		outputDir = outputDir
		)

	sim = wholecell.sim.simulation.Simulation(**simOpts)

	sim.run()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("submissionTime", help = "Time of job submission", type = str)

	args = parser.parse_args().__dict__

	main(args["submissionTime"])

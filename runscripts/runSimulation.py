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

import models.ecoli.sim.sim_definition
import wholecell.sim.simulation
import argparse
import os
import json

def main(submissionTime):

	simOpts = models.ecoli.sim.sim_definition.getSimOptsFromEnvVars(
		["outputDir", "logToDisk", "overwriteExistingFiles"]
		)

	# Define simulation options specific to this script
	simOpts["outputDir"] = os.path.join(
		"out",
		"%s" % submissionTime,
		"%06d" % simOpts["seed"],
		"simOut",
		)
	simOpts["logToDisk"] = True
	simOpts["overwriteExistingFiles"] = False

	# Check that we're setting all arguments (in case more have been added, etc)
	assert (
		set(simOpts.keys()) ==
		set(models.ecoli.sim.sim_definition.SIM_KWARG_DEFAULTS.keys())
		), "Need to set all keyword arguments in runSimulation.py"

	sim = wholecell.sim.simulation.Simulation(models.ecoli.sim.sim_definition.SimDefinition(**simOpts))

	sim.run()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("submissionTime", help = "Time of job submission", type = str)

	args = parser.parse_args().__dict__

	main(args["submissionTime"])

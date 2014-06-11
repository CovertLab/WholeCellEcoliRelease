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

import wholecell.sim.sim_definition
import wholecell.sim.simulation
import argparse
import os
import json

def main(submissionTime):

	# We use this to check if any undefined WC_* environmental variables
	# were accidentally specified by the user
	wcEnvVars = [x for x in os.environ if x.startswith("WC_")]

	optionsAndEnvVars = dict(
		seed = ("WC_SEED", int),
		states = ("WC_STATES", json.loads),
		processes = ("WC_PROCESSES", json.loads),
		listeners = ("WC_LISTENERS", json.loads),
		hooks = ("WC_HOOKS", json.loads),
		lengthSec = ("WC_LENGTHSEC", int),
		timeStepSec = ("WC_TIMESTEPSEC", float),
		logToShell = ("WC_LOGTOSHELL", json.loads),
		shellColumnHeaders = ("WC_SHELLCOLUMNSHEADERS", json.loads),
		logToDiskEvery = ("WC_LOGTODISKEVERY", int),
		kbLocation = ("WC_KBLOCATION", json.loads)
		)

	simOpts = {}

	# Get simulation options from environmental variables
	for option, (envVar, handler) in optionsAndEnvVars.iteritems():
		if os.environ.has_key(envVar) and len(os.environ[envVar]):
			simOpts[option] = handler(os.environ[envVar])
			wcEnvVars.remove(envVar)
		else:
			simOpts[option] = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS[option]

	# Check for extraneous environmental variables (probably typos by the user)
	assert (len(wcEnvVars) == 0), (
		"The following WC_* environmental variables were specified but " +
		"have no defined function: %s" % wcEnvVars
		)

	# Define simulation options specific to this script
	simOpts["outputDir"] = os.path.join(
		"out",
		"simOut",
		"%s" % submissionTime,
		"%06d" % simOpts["seed"]
		)
	simOpts["logToDisk"] = True
	simOpts["overwriteExistingFiles"] = False

	# Check that we're setting all arguments (in case more have been added, etc)
	assert (
		set(simOpts.keys()) ==
		set(wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS.keys())
		), "Need to set all keyword arguments in runSimulation.py"

	sim = wholecell.sim.simulation.Simulation(**simOpts)

	sim.run()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("submissionTime", help = "Time of job submission", type = str)

	args = parser.parse_args().__dict__

	main(args["submissionTime"])

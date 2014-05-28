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

	seed = 0
	if os.environ.has_key("ARRAY_ID"):
		seed = int(os.environ["ARRAY_ID"]) - 1

	states = wholecell.sim.sim_definition.DEFAULT_STATES
	if os.environ.has_key("WC_STATES") and len(os.environ["WC_STATES"]):
		print os.environ["WC_STATES"]
		states = json.loads(os.environ["WC_STATES"])

	processes = wholecell.sim.sim_definition.DEFAULT_PROCESSES
	if os.environ.has_key("WC_PROCESSES") and len(os.environ["WC_PROCESSES"]):
		processes = json.loads(os.environ["WC_PROCESSES"])

	listeners = wholecell.sim.sim_definition.DEFAULT_LISTENERS
	if os.environ.has_key("WC_LISTENERS") and len(os.environ["WC_LISTENERS"]):
		listeners = json.loads(os.environ["WC_LISTENERS"])

	hooks = wholecell.sim.sim_definition.DEFAULT_HOOKS
	if os.environ.has_key("WC_HOOKS") and len(os.environ["WC_HOOKS"]):
		hooks = json.loads(os.environ["WC_HOOKS"])

	lengthSec = wholecell.sim.sim_definition.DEFAULT_LENGTH
	if os.environ.has_key("WC_LENGTHSEC") and len(os.environ["WC_LENGTHSEC"]):
		lengthSec = int(os.environ["WC_LENGTHSEC"])

	timeStepSec = wholecell.sim.sim_definition.DEFAULT_TIME_STEP
	if os.environ.has_key("WC_TIMESTEPSEC") and len(os.environ["WC_TIMESTEPSEC"]):
		timeStepSec = float(timeStepSec)

	logToShell = True
	if os.environ.has_key("WC_LOGTOSHELL") and len(os.environ["WC_LOGTOSHELL"]):
		logToShell = json.loads(os.environ["WC_LOGTOSHELL"])

	logToDiskEvery = None
	if os.environ.has_key("WC_LOGTODISKEVERY") and len(os.environ["WC_LOGTODISKEVERY"]):
		logToDiskEvery = int(os.environ["WC_LOGTODISKEVERY"])

	rebuildKB = True
	if os.environ.has_key("WC_REBUILDKB") and len(os.environ["WC_REBUILDKB"]):
		rebuildKB = json.loads(os.environ["WC_REBUILDKB"])

	outputDir = os.path.join("out", "simOut", "%s" % submissionTime, "%06d" % seed)

	simOpts = dict(
		states = states,
		processes = processes,
		listeners = listeners,
		hooks = hooks,
		lengthSec = lengthSec,
		timeStepSec = timeStepSec,
		seed = seed,
		logToShell = logToShell,
		logToDisk = True,
		logToDiskEvery = logToDiskEvery,
		outputDir = outputDir
		)

	sim = wholecell.sim.simulation.Simulation(**simOpts)

	sim.run()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("submissionTime", help = "Time of job submission", type = str)

	args = parser.parse_args().__dict__

	main(args["submissionTime"])
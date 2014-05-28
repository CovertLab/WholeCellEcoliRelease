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
import os
import json

def main():

	# Run simulation with options specified in a json file
	if os.environ.has_key("WC_FILE") and len(os.environ["WC_FILE"]):
		sim = wholecell.sim.simulation.Simulation.initFromFile(
			os.environ["WC_FILE"]
			)
		sim.run()
		return

	# Get simulation options from environmental variables
	seed = 0
	if os.environ.has_key("WC_SEED") and len(os.environ["WC_SEED"]):
		seed = int(os.environ["WC_SEED"])

	states = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["states"]
	if os.environ.has_key("WC_STATES") and len(os.environ["WC_STATES"]):
		states = json.loads(os.environ["WC_STATES"])

	processes = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["processes"]
	if os.environ.has_key("WC_PROCESSES") and len(os.environ["WC_PROCESSES"]):
		processes = json.loads(os.environ["WC_PROCESSES"])

	listeners = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["listeners"]
	if os.environ.has_key("WC_LISTENERS") and len(os.environ["WC_LISTENERS"]):
		listeners = json.loads(os.environ["WC_LISTENERS"])

	hooks = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["hooks"]
	if os.environ.has_key("WC_HOOKS") and len(os.environ["WC_HOOKS"]):
		hooks = json.loads(os.environ["WC_HOOKS"])

	lengthSec = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["lengthSec"]
	if os.environ.has_key("WC_LENGTHSEC") and len(os.environ["WC_LENGTHSEC"]):
		lengthSec = int(os.environ["WC_LENGTHSEC"])

	timeStepSec = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["timeStepSec"]
	if os.environ.has_key("WC_TIMESTEPSEC") and len(os.environ["WC_TIMESTEPSEC"]):
		timeStepSec = float(timeStepSec)

	logToShell = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["logToShell"]
	if os.environ.has_key("WC_LOGTOSHELL") and len(os.environ["WC_LOGTOSHELL"]):
		logToShell = json.loads(os.environ["WC_LOGTOSHELL"])

	logToDisk = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["logToDisk"]
	if os.environ.has_key("WC_LOGTODISK") and len(os.environ["WC_LOGTODISK"]):
		logToDisk = json.loads(os.environ["WC_LOGTODISK"])

	outputDir = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["outputDir"]
	if os.environ.has_key("WC_OUTPUTDIR") and len(os.environ["WC_OUTPUTDIR"]):
		outputDir = json.loads(os.environ["WC_OUTPUTDIR"])

	overwriteExistingFiles = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["overwriteExistingFiles"]
	if os.environ.has_key("WC_OVERWRITEEXISTINGFILES") and len(os.environ["WC_OVERWRITEEXISTINGFILES"]):
		overwriteExistingFiles = json.loads(os.environ["WC_OVERWRITEEXISTINGFILES"])

	logToDiskEvery = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["logToDiskEvery"]
	if os.environ.has_key("WC_LOGTODISKEVERY") and len(os.environ["WC_LOGTODISKEVERY"]):
		logToDiskEvery = int(os.environ["WC_LOGTODISKEVERY"])

	rebuildKB = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["rebuildKB"]
	if os.environ.has_key("WC_REBUILDKB") and len(os.environ["WC_REBUILDKB"]):
		rebuildKB = json.loads(os.environ["WC_REBUILDKB"])

	simOpts = dict(
		states = states,
		processes = processes,
		listeners = listeners,
		hooks = hooks,
		lengthSec = lengthSec,
		timeStepSec = timeStepSec,
		seed = seed,
		logToShell = logToShell,
		logToDisk = logToDisk,
		outputDir = outputDir,
		overwriteExistingFiles = overwriteExistingFiles,
		logToDiskEvery = logToDiskEvery,
		rebuildKB = rebuildKB
		)

	assert (
		set(simOpts.keys()) ==
		set(wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS.keys())
		), "Need to set all keyword arguments in runSimulation.py"

	sim = wholecell.sim.simulation.Simulation(**simOpts)

	sim.run()

if __name__ == '__main__':
	main()

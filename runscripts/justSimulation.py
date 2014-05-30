#!/usr/bin/env python

"""
justSimulation.py

Runs a simulation.  Called from justSimulation.sh

"""

from __future__ import division

import wholecell.sim.sim_definition
import wholecell.sim.simulation
import os
import json

def main():

	# Get simulation options from environmental variables
	seed = 0
	if os.environ.has_key("WC_SEED") and len(os.environ["WC_SEED"]):
		seed = int(os.environ["WC_SEED"])

	# We use this to check if any undefined WC_* environmental variables
	# were accidentally specified by the user
	wcEnvVars = [x for x in os.environ if x.startswith("WC_")]

	states = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["states"]
	if os.environ.has_key("WC_STATES") and len(os.environ["WC_STATES"]):
		states = json.loads(os.environ["WC_STATES"])
		wcEnvVars.remove("WC_STATES")

	processes = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["processes"]
	if os.environ.has_key("WC_PROCESSES") and len(os.environ["WC_PROCESSES"]):
		processes = json.loads(os.environ["WC_PROCESSES"])
		wcEnvVars.remove("WC_PROCESSES")

	listeners = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["listeners"]
	if os.environ.has_key("WC_LISTENERS") and len(os.environ["WC_LISTENERS"]):
		listeners = json.loads(os.environ["WC_LISTENERS"])
		wcEnvVars.remove("WC_LISTENERS")

	hooks = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["hooks"]
	if os.environ.has_key("WC_HOOKS") and len(os.environ["WC_HOOKS"]):
		hooks = json.loads(os.environ["WC_HOOKS"])
		wcEnvVars.remove("WC_HOOKS")

	lengthSec = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["lengthSec"]
	if os.environ.has_key("WC_LENGTHSEC") and len(os.environ["WC_LENGTHSEC"]):
		lengthSec = int(os.environ["WC_LENGTHSEC"])
		wcEnvVars.remove("WC_LENGTHSEC")

	timeStepSec = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["timeStepSec"]
	if os.environ.has_key("WC_TIMESTEPSEC") and len(os.environ["WC_TIMESTEPSEC"]):
		timeStepSec = float(timeStepSec)
		wcEnvVars.remove("WC_TIMESTEPSEC")

	logToShell = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["logToShell"]
	if os.environ.has_key("WC_LOGTOSHELL") and len(os.environ["WC_LOGTOSHELL"]):
		logToShell = json.loads(os.environ["WC_LOGTOSHELL"])
		wcEnvVars.remove("WC_LOGTOSHELL")

	logToDisk = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["logToDisk"]
	if os.environ.has_key("WC_LOGTODISK") and len(os.environ["WC_LOGTODISK"]):
		logToDisk = int(os.environ["WC_LOGTODISK"])
		wcEnvVars.remove("WC_LOGTODISK")

	outputDir = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["outputDir"]
	if os.environ.has_key("WC_OUTPUTDIR") and len(os.environ["WC_OUTPUTDIR"]):
		outputDir = int(os.environ["WC_OUTPUTDIR"])
		wcEnvVars.remove("WC_OUTPUTDIR")

	overwriteExistingFiles = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["overwriteExistingFiles"]
	if os.environ.has_key("WC_OVERWRITEEXISTINGFILES") and len(os.environ["WC_OVERWRITEEXISTINGFILES"]):
		overwriteExistingFiles = int(os.environ["WC_OVERWRITEEXISTINGFILES"])
		wcEnvVars.remove("WC_OVERWRITEEXISTINGFILES")

	logToDiskEvery = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["logToDiskEvery"]
	if os.environ.has_key("WC_LOGTODISKEVERY") and len(os.environ["WC_LOGTODISKEVERY"]):
		logToDiskEvery = int(os.environ["WC_LOGTODISKEVERY"])
		wcEnvVars.remove("WC_LOGTODISKEVERY")

	kbLocation = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["kbLocation"]
	if os.environ.has_key("WC_KBLOCATION") and len(os.environ["WC_KBLOCATION"]):
		kbLocation = json.loads(os.environ["WC_KBLOCATION"])
		wcEnvVars.remove("WC_KBLOCATION")

	assert (len(wcEnvVars) == 0), (
		"The following WC_* environmental variables were specified but " +
		"have no defined function: %s" % wcEnvVars
		)

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
		kbLocation = kbLocation
		)

	assert (
		set(simOpts.keys()) ==
		set(wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS.keys())
		), "Need to set all keyword arguments in runSimulation.py"

	sim = wholecell.sim.simulation.Simulation(**simOpts)

	sim.run()

if __name__ == '__main__':
	main()

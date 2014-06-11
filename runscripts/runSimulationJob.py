#!/usr/bin/env python

"""
runSimulationJob.py

Runs a simulation, called from a job script.

Options are passed using environmental variables.

"""

from __future__ import division

import wholecell.sim.sim_definition
import wholecell.sim.simulation
import argparse
import os
import json

def main(submissionTime):

	seed = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["seed"]
	if os.environ.has_key("ARRAY_ID"):
		seed = int(os.environ["ARRAY_ID"]) - 1

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

	shellColumnHeaders = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["shellColumnHeaders"]
	if os.environ.has_key("WC_SHELLCOLUMNSHEADERS") and len(os.environ["WC_SHELLCOLUMNSHEADERS"]):
		shellColumnHeaders = json.loads(os.environ["WC_SHELLCOLUMNSHEADERS"])
		wcEnvVars.remove("WC_SHELLCOLUMNSHEADERS")

	logToDiskEvery = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["logToDiskEvery"]
	if os.environ.has_key("WC_LOGTODISKEVERY") and len(os.environ["WC_LOGTODISKEVERY"]):
		logToDiskEvery = int(os.environ["WC_LOGTODISKEVERY"])
		wcEnvVars.remove("WC_LOGTODISKEVERY")

	rebuildKB = wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS["rebuildKB"]
	if os.environ.has_key("WC_REBUILDKB") and len(os.environ["WC_REBUILDKB"]):
		rebuildKB = json.loads(os.environ["WC_REBUILDKB"])
		wcEnvVars.remove("WC_REBUILDKB")

	assert (len(wcEnvVars) == 0), (
		"The following WC_* environmental variables were specified but " +
		"have no defined function: %s" % wcEnvVars
		)

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
		shellColumnHeaders = shellColumnHeaders,
		logToDisk = True,
		outputDir = outputDir,
		overwriteExistingFiles = False,
		logToDiskEvery = logToDiskEvery,
		rebuildKB = rebuildKB
		)

	assert (
		set(simOpts.keys()) ==
		set(wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS.keys())
		), "Need to set all keyword arguments in runSimulationJob.py"

	sim = wholecell.sim.simulation.Simulation(**simOpts)

	sim.run()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("submissionTime", help = "Time of job submission", type = str)

	args = parser.parse_args().__dict__

	main(args["submissionTime"])
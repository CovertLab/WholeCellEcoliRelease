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

	simOpts = wholecell.sim.sim_definition.getSimOptsFromEnvVars()

	# Check that we're setting all arguments (in case more have been added, etc)
	assert (
		set(simOpts.keys()) ==
		set(wholecell.sim.sim_definition.SIM_KWARG_DEFAULTS.keys())
		), "Need to set all keyword arguments in runSimulation.py"

	sim = wholecell.sim.simulation.Simulation(**simOpts)

	sim.run()

if __name__ == '__main__':
	main()

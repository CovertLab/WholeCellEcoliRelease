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
import sys

DEFAULT_SIM = dict(
	seed = 10,
	lengthSec = 10,
	includedProcesses = ['ToyReplication']
	)

def main():
	# TODO: argument parsing

	nArgs = len(sys.argv)

	if nArgs == 1:
		# Use default parameters
		sim = wholecell.sim.simulation.Simulation(**DEFAULT_SIM)

	elif nArgs == 2:
		# Attempt to parse from a json file
		sim = wholecell.sim.simulation.Simulation.initFromFile(sys.argv[1])

	sim.run()

if __name__ == '__main__':
	main()

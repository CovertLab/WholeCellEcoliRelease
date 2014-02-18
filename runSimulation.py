#!/usr/bin/env python

"""
runSimulation
Runs and logs whole-cell simulation

Example:
>>> from runSimulation import *
>>> runSimulation()

Example:
~/parWholeCellPy$ python runSimulation.py

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/10/2013
"""

import wholecell.sim.simulation
import sys

DEFAULT_SIM = dict(
	seed = 10,
	lengthSec = 100,
	logToDisk = True
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

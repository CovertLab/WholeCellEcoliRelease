#!/usr/bin/env python

"""
loadSimulation.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/18/2013
"""

import wholecell.sim.simulation
import sys

def main():
	# TODO: add optional arguments

	sim = wholecell.sim.simulation.Simulation.loadSimulation(
		sys.argv[1], int(sys.argv[2]),
		)

	sim.run()

if __name__ == '__main__':
	main()

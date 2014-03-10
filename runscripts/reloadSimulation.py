#!/usr/bin/env python

"""
reloadSimulation.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/18/2013
"""

from __future__ import division

import wholecell.sim.simulation
import sys

def main():
	# TODO: add optional arguments
	# Takes in two arguments from a command line call:
	# Argument 1: Directory containing the saved simulation files
	# Argument 2: The time step to load from
	#
	# Example call: python reloadSimulation.py /out/savedsim 0
	
	sim = wholecell.sim.simulation.Simulation.loadSimulation(
		sys.argv[1], int(sys.argv[2]),
		)

	sim.run()

if __name__ == '__main__':
	main()

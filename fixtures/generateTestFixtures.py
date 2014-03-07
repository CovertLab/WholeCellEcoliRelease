#!/usr/bin/env python

"""
generateTestFixtures
Generates fixtures for whole-cell tests

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/5/2013
"""

import os
import cPickle

import wholecell.reconstruction.knowledgebase
import wholecell.sim.simulation
import wholecell.utils.knowledgebase_fixture_manager

import wholecell.utils.config_test

def main():
	# Create output directory
	outDir = wholecell.utils.config.TEST_FIXTURE_DIR
	if not os.path.exists(outDir):
		os.makedirs(outDir)

	# Construct KB
	wholecell.utils.knowledgebase_fixture_manager.cacheKnowledgeBase(outDir)

	# Construct simulation
	sim = wholecell.sim.simulation.Simulation(
		seed = 1
		)

	cPickle.dump(sim, open(os.path.join(outDir, "Simulation.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)
	# TODO: save/load a simulation's initial step instead of pickling the simulation

if __name__ == '__main__':
	main()

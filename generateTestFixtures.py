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

def main():
	# Create output directory
	outDir = "data/fixtures"
	if not os.path.exists(outDir):
		os.makedirs(outDir)

	# Construct KB
	kb = wholecell.reconstruction.knowledgebase.KnowledgeBase(dataFileDir = "data/parsed/", seqFileName = "data/raw/sequence.txt")
	cPickle.dump(kb, open(os.path.join(outDir, "KnowledgeBase.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)

	# Construct simulation
	sim = wholecell.sim.simulation.Simulation(
		seed = 1
		)

	cPickle.dump(sim, open(os.path.join(outDir, "Simulation.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
	main()

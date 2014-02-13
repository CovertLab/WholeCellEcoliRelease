#!/usr/bin/env python

"""
generateTestFixtures
Generates fixtures for whole-cell tests

Example:
>>> from generateTestFixtures import *
>>> generateTestFixtures()

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/5/2013
"""

import os
import cPickle

# We need the following two methods so that we can pickle instancemethods
# Also need the copy_reg.pickle() line in generateTestFixtures()
# Borrowed from: http://mail.python.org/pipermail/python-list/2006-October/367078.html
# (Thanks Steven Bethard!)

def generateTestFixtures():
	import wholecell.kb.KnowledgeBase
	import wholecell.sim.Simulation
	import wholecell.util.Fitter

	# Create output directory
	outDir = "data/fixtures"
	if not os.path.exists(outDir):
		os.makedirs(outDir)

	# Construct KB
	kb = wholecell.kb.KnowledgeBase.KnowledgeBase(dataFileDir = "data/parsed/", seqFileName = "data/raw/sequence.txt")
	cPickle.dump(kb, open(os.path.join(outDir, "KnowledgeBase.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)

	# Construct simulation
	sim = wholecell.sim.Simulation.Simulation()
	sim.initialize(kb)
	sim.setOptions({"seed": 1})
	wholecell.util.Fitter.Fitter.FitSimulation(sim, kb)
	cPickle.dump(sim, open(os.path.join(outDir, "Simulation.cPickle"), "wb"), protocol = cPickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
	generateTestFixtures()

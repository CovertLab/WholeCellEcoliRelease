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

import os
import cPickle

import wholecell.loggers.disk
import wholecell.loggers.shell
import wholecell.sim.simulation
import wholecell.reconstruction.knowledgebase
import wholecell.utils.fitter

import weakref

KB_PATH = os.path.join('data', 'fixtures', 'KnowledgeBase.cPickle')

def runSimulation(reconstructKB = False, fitSimulation = True,
		useShellLogger = True, useDiskLogger = False, outDir = None,
		simOpts = None):

	# Instantiate knowledge base
	if reconstructKB or not os.path.exists(KB_PATH):
		kb = wholecell.reconstruction.knowledgebase.KnowledgeBase(
			dataFileDir = "data/parsed", seqFileName = "data/raw/sequence.txt"
			)

		cPickle.dump(kb, open(KB_PATH, "wb"),
			protocol = cPickle.HIGHEST_PROTOCOL)

	else:
		kb = cPickle.load(open(KB_PATH, "rb"))

	kbWeakRef = weakref.ref(kb)

	# Set up simulation
	sim = wholecell.sim.simulation.Simulation()

	sim.initialize(kb)

	if simOpts:
		sim.setOptions(simOpts)

	if fitSimulation:
		wholecell.utils.fitter.Fitter.FitSimulation(sim, kb)

	del kb

	assert kbWeakRef() is None, 'Failed to release knowledge base from memory.'

	# Instantiate loggers
	if useShellLogger:
		sim.loggerAdd(wholecell.loggers.shell.Shell())

	if useDiskLogger:
		sim.loggerAdd(
			wholecell.loggers.disk.Disk(outDir = outDir)
			)

	# Run simulation
	sim.run()

if __name__ == '__main__':
	runSimulation(
		simOpts = {
			'seed':10,
			'lengthSec':100
			},
		fitSimulation = True,
		useDiskLogger = False,
		)

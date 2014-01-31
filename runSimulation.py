#!/usr/bin/env python

"""
runSimulation
Runs and logs whole-cell simulation

Example:
>>> from runSimulation import *
>>> runSimulation()

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/10/2013
"""

import os
import cPickle

import numpy

import wholecell.sim.logger.Disk
import wholecell.sim.logger.Shell
import wholecell.sim.Simulation
import wholecell.kb.KnowledgeBase
import wholecell.util.Fitter

import weakref

KB_PATH = os.path.join('data', 'fixtures', 'KnowledgeBase.cPickle')

DEFAULT_OPTIONS = {
	'seed':10,
	'lengthSec':100
	}

def runSimulation(reconstructKB = False, fitSimulation = True,
		useShellLogger = True, useDiskLogger = False, outFile = None,
		simOpts = None):

	# Instantiate knowledge base
	if reconstructKB or not os.path.exists(KB_PATH):
		kb = wholecell.kb.KnowledgeBase.KnowledgeBase(
			dataFileDir = "data/parsed", seqFileName = "data/raw/sequence.txt"
			)

		cPickle.dump(kb, open(KB_PATH, "wb"),
			protocol = cPickle.HIGHEST_PROTOCOL)

	else:
		kb = cPickle.load(open(KB_PATH, "rb"))

	kbWeakRef = weakref.ref(kb)

	# Set up simulation
	sim = wholecell.sim.Simulation.Simulation(kb)
	if simOpts:
		sim.setOptions(simOpts)

	if fitSimulation:
		wholecell.util.Fitter.Fitter.FitSimulation(sim, kb)

	del kb

	assert kbWeakRef() is None, 'Failed to release knowledge base from memory.'

	# Instantiate loggers
	loggers = []

	if useShellLogger:
		loggers.append(wholecell.sim.logger.Shell.Shell())

	if useDiskLogger:
		if outFile is None:
			raise Exception('No output directory provided.')

		loggers.append(
			wholecell.sim.logger.Disk.Disk(outFile = outFile, overwrite = True)
			)

	# Run simulation
	sim.run(loggers)

if __name__ == '__main__':
	runSimulation(simOpts = DEFAULT_OPTIONS, useDiskLogger = True, outFile = 'out/stateFile.hdf')

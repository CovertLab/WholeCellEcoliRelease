
from __future__ import division

import cPickle
import os
import json
import pkg_resources

import numpy as np

import wholecell.sim.simulation as wcSimulation
import wholecell.loggers.disk as wcDisk

from mpi4py import MPI

COMM = MPI.COMM_WORLD

import wholecell.utils.config
TEST_FIXTURE_DIR = wholecell.utils.config.TEST_FIXTURE_DIR
KB_PATH = os.path.join(TEST_FIXTURE_DIR,'KnowledgeBase.cPickle')

def main():
	fixtureOpts = json.load(open('fixture_tmp.json'))

	# TODO: make splitting the option file more modular?
	simOpts = fixtureOpts.copy()
	simOpts.pop('nSeeds')
	simOpts['logToDisk'] = True
	simOpts['overwriteExistingFiles'] = True
	simOpts['autoRun'] = True

	N_SEEDS = fixtureOpts['nSeeds']
	N_PROC = COMM.size # the number of processes
	SEND_COUNTS = (
		(N_SEEDS // N_PROC) * np.ones(N_PROC, dtype = int)
		+ (np.arange(N_PROC, dtype =int) < N_SEEDS % N_PROC)
		) # the number of seeds each process receives
	DISPLACEMENTS = np.hstack([
		np.zeros(1),
		np.cumsum(SEND_COUNTS)[:-1]
		]) # the starting position for each process's seeds

	##################################################
	# MPI seed scatter logic
	##################################################
	mySeeds = np.zeros(SEND_COUNTS[COMM.rank]) # buffer for the seeds received
	allSeeds = np.arange(N_SEEDS, dtype = float) # the list of seeds

	COMM.Scatterv(
		[allSeeds, tuple(SEND_COUNTS), tuple(DISPLACEMENTS), MPI.DOUBLE],
		mySeeds
		) # scatter data from the root (by default) to each child process

	print "Rank %d mySeeds: %s" % (COMM.rank, str(mySeeds))
	##################################################

	for seed in mySeeds.astype('int'):
		simOpts['outputDir'] = os.path.join(
			fixtureOpts['outputDir'], 'sim{}'.format(seed)
			)

		print 'Running simulation #{}'.format(seed)
		sim = wcSimulation.Simulation(**simOpts)


if __name__ == '__main__':
	main()

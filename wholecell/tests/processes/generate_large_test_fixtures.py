import cPickle
import os
import json
import pkg_resources

import numpy

import wholecell.utils.fitter as wcFitter
import wholecell.sim.simulation as wcSimulation
import wholecell.loggers.disk as wcDisk

from mpi4py import MPI

COMM = MPI.COMM_WORLD

KB_PATH = pkg_resources.resource_filename('data','fixtures/KnowledgeBase.cPickle')

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
		(N_SEEDS // N_PROC) * numpy.ones(N_PROC, dtype = int)
		+ (numpy.arange(N_PROC, dtype =int) < N_SEEDS % N_PROC)
		) # the number of seeds each process receives
	DISPLACEMENTS = numpy.hstack([
		numpy.zeros(1),
		numpy.cumsum(SEND_COUNTS)[:-1]
		]) # the starting position for each process's seeds

	##################################################
	# MPI seed scatter logic
	##################################################
	mySeeds = numpy.zeros(SEND_COUNTS[COMM.rank]) # buffer for the seeds received
	allSeeds = numpy.arange(N_SEEDS, dtype = float) # the list of seeds

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

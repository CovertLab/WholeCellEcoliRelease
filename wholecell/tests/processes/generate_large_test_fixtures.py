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

def runSimulations():
	infile = open('fixture_tmp.json', 'r')
	fixture_opts = json.loads(infile.readline())
	infile.close()

	LENGTH_SEC = fixture_opts['length_sec']
	KB = cPickle.load(open(KB_PATH, "rb"))

	N_SEEDS = fixture_opts['n_seeds']
	N_PROC = COMM.size # the number of processes
	SEND_COUNTS = (
		(N_SEEDS // N_PROC) * numpy.ones(N_PROC, dtype = int)
		+ (numpy.arange(N_PROC, dtype =int) < N_SEEDS % N_PROC)
		) # the number of seeds each process receives
	DISPLACEMENTS = numpy.hstack([
		numpy.zeros(1),
		numpy.cumsum(SEND_COUNTS)[:-1]
		]) # the starting position for each process's seeds

	logDir = fixture_opts['fixture_dir']

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
		sim = wcSimulation.Simulation(fixture_opts['processes'], fixture_opts['free_molecules'])

		sim.initialize(KB)

		wcFitter.Fitter.FitSimulation(sim, KB)

		sim.setOptions({"lengthSec":LENGTH_SEC, "seed":seed})

		sim.loggerAdd(
			wcDisk.Disk(
				os.path.join(logDir, 'sim{}'.format(seed)),
				True
				)
			)

		print 'Running simulation #{}'.format(seed)
		sim.run()


def main():
	runSimulations()

if __name__ == '__main__':
	main()
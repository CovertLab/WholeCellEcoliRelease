import cPickle
import os

import numpy

import wholecell.util.Fitter as wcFitter
import wholecell.sim.Simulation as wcSimulation
import wholecell.sim.logger.Disk as wcDisk

from mpi4py import MPI

COMM = MPI.COMM_WORLD

# TODO: establish naming scheme for test simulation output
# TODO: cache fit simulation prior to running/logging
# TODO: save fit parameters in a cached knowledge base object

KB_PATH = os.path.join('data', 'fixtures', 'KnowledgeBase.cPickle')

LENGTH_SEC = 500
KB = cPickle.load(open(KB_PATH, "rb"))

N_SEEDS = 8
N_PROC = COMM.size # the number of processes
SEND_COUNTS = (
	(N_SEEDS // N_PROC) * numpy.ones(N_PROC, dtype = int)
	+ (numpy.arange(N_PROC, dtype =int) < N_SEEDS % N_PROC)
	) # the number of seeds each process receives
DISPLACEMENTS = numpy.hstack([
	numpy.zeros(1),
	numpy.cumsum(SEND_COUNTS)[:-1]
	]) # the starting position for each process's seeds

def runSimulations(testDir, processes, freeMolecules):
	logDir = os.path.join('out', 'test', testDir)

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
		sim = wcSimulation.Simulation(processes, freeMolecules)

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
	ntpCounts = 1e6
	initEnzCnts = 2000.

	# Tests for Transcription-only simulations
	runSimulations(
		testDir = 'Test_Transcription',
		processes = ["Transcription"],
		freeMolecules = [
			["ATP[c]", ntpCounts],
			["UTP[c]", ntpCounts],
			["CTP[c]", ntpCounts],
			["GTP[c]", ntpCounts],
			["EG10893-MONOMER[c]", initEnzCnts],
			["RPOB-MONOMER[c]", initEnzCnts],
			["RPOC-MONOMER[c]", initEnzCnts],
			["RPOD-MONOMER[c]", initEnzCnts]
		]
		)

if __name__ == '__main__':
	main()
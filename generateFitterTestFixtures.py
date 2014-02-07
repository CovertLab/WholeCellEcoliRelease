'''

Let's replace this file with something better - John

'''

import cPickle
import os

import numpy

import wholecell.util.Fitter as wcFitter
import wholecell.sim.Simulation as wcSimulation
import wholecell.sim.logger.Disk as wcDisk

from mpi4py import MPI

comm = MPI.COMM_WORLD

# TODO: make parallel
# TODO: make fixture generation more modular
# TODO: establish naming scheme for test simulation output
# TODO: cache fit simulation prior to running/logging
# TODO: save fit parameters in a cached knowledge base object

KB_PATH = os.path.join('data', 'fixtures', 'KnowledgeBase.cPickle')

ntpCounts = 1e6
initEnzCnts = 2000.
# initRnaCnts = 0.
# T_d = 3600.
lengthSec = 500
nSeeds = 8

logDir = os.path.join('out', 'test', 'rnaProduction')

kb = cPickle.load(open(KB_PATH, "rb"))

##################################################
# MPI seed scatter logic
##################################################

nProc = comm.size # the number of processes
sendcounts = (
	(nSeeds // nProc) * numpy.ones(nProc, dtype = int)
	+ (numpy.arange(nProc, dtype =int) < nSeeds % nProc)
	) # the number of seeds each process receives
displacements = numpy.hstack([
	numpy.zeros(1),
	numpy.cumsum(sendcounts)[:-1]
	]) # the starting position for each process's seeds
mySeeds = numpy.zeros(sendcounts[comm.rank]) # buffer for the seeds received
allSeeds = numpy.arange(nSeeds, dtype = float) # the list of seeds

comm.Scatterv(
	[allSeeds, tuple(sendcounts), tuple(displacements), MPI.DOUBLE],
	mySeeds
	) # scatter data from the root (by default) to each child process

print "Rank %d mySeeds: %s" % (comm.rank, str(mySeeds))
##################################################

for seed in mySeeds.astype('int'):
	sim = wcSimulation.Simulation(
		["Transcription"],
		[
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

	sim.initialize(kb)

	wcFitter.Fitter.FitSimulation(sim, kb)

	sim.setOptions({"lengthSec":lengthSec, "seed":seed})

	sim.loggerAdd(
		wcDisk.Disk(
			os.path.join(logDir, 'sim{}'.format(seed)),
			True
			)
		)

	print 'Running simulation #{}'.format(seed)
	sim.run()

import cPickle
import os

import numpy as np

import wholecell.utils.fitter as wcFitter
import wholecell.sim.simulation as wcSimulation
import wholecell.loggers.disk as wcDisk

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
	(N_SEEDS // N_PROC) * np.ones(N_PROC, dtype = int)
	+ (np.arange(N_PROC, dtype =int) < N_SEEDS % N_PROC)
	) # the number of seeds each process receives
DISPLACEMENTS = np.hstack([
	np.zeros(1),
	np.cumsum(SEND_COUNTS)[:-1]
	]) # the starting position for each process's seeds

def runSimulations(testDir, processes, freeMolecules):
	logDir = os.path.join('out', 'test', testDir)

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
	# Tests for Transcription-only simulations
	ntpCounts = 1e6
	initEnzCnts = 2000.
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

	# Tests for Transcription + RnaDegradation/Maturation-only simulations
	initRnapCnts = 1000.
	initRnaseCnts = 1000.

	ntpCounts = 1e6
	h2oCounts = 1e6

	# TODO: increase simulation time to 3600
	# NOTE: test only runs and plots one simulation

	runSimulations(
		testDir = 'Test_Transcription_RnaDegradation',
		processes = ["Transcription", "RnaDegradation", "RnaMaturation"],
		freeMolecules = [
			["ATP[c]", ntpCounts],
			["UTP[c]", ntpCounts],
			["CTP[c]", ntpCounts],
			["GTP[c]", ntpCounts],
			["H2O[c]", h2oCounts],
			["EG10893-MONOMER[c]", initRnapCnts],
			["RPOB-MONOMER[c]", initRnapCnts],
			["RPOC-MONOMER[c]", initRnapCnts],
			["RPOD-MONOMER[c]", initRnapCnts],
			["EG11259-MONOMER[c]", initRnaseCnts]
		]
		)

	# Tests for Transcription + RnaDegradation/Maturation + Translation + ProteinMaturation-only simulations
	initRnapCnts = 1000.
	initRnaseCnts = 1000.
	initRibCnts = 1000. / 7 # seven isozymes?

	ntpCounts = 1e6
	h2oCounts = 1e6
	aaCounts = 1e8

	# TODO: increase simulation time to 3600

	from wholecell.processes.translation import enzIDs as translationEnzymes
	from wholecell.processes.translation import aaIDs

	runSimulations(
		testDir = 'Test_Transcription_RnaDegradation_Translation',
		processes = ["Transcription", "RnaDegradation", "RnaMaturation", 'Translation', 'ProteinMaturation'],
		freeMolecules = [
			["ATP[c]", ntpCounts],
			["UTP[c]", ntpCounts],
			["CTP[c]", ntpCounts],
			["GTP[c]", ntpCounts],
			["H2O[c]", h2oCounts],
			["EG10893-MONOMER[c]", initRnapCnts],
			["RPOB-MONOMER[c]", initRnapCnts],
			["RPOC-MONOMER[c]", initRnapCnts],
			["RPOD-MONOMER[c]", initRnapCnts],
			["EG11259-MONOMER[c]", initRnaseCnts]
		] + [
			[enzID, initRibCnts] for enzID in translationEnzymes
		] + [
			[aaID, aaCounts] for aaID in aaIDs
		]
		)

	# Tests for RnaDegradation-only simulations
	initRnaseCnts = 1000.
	# NOTE: in the original test, this set every RNA counts to 10000
	# TODO: implement the ability to set initial counts, or remove this note

	h2oCounts = 1e6

	from wholecell.processes.Translation import enzIDs as translationEnzymes
	from wholecell.processes.Translation import aaIDs

	runSimulations(
		testDir = 'Test_RnaDegradation',
		processes = ["RnaDegradation"],
		freeMolecules = [
			["H2O[c]", h2oCounts],
			["EG11259-MONOMER[c]", initRnaseCnts],
		]
		)

	# Tests for Translation-only simulations
	initRibCnts = 1000. / 7 # seven isozymes?

	h2oCounts = 1e6
	aaCounts = 1e8

	# TODO: provide with mRNAs

	from wholecell.processes.Translation import enzIDs as translationEnzymes
	from wholecell.processes.Translation import aaIDs

	runSimulations(
		testDir = 'Test_Translation',
		processes = ['Translation'],
		freeMolecules = [
			[enzID, initRibCnts] for enzID in translationEnzymes
		] + [
			[aaID, aaCounts] for aaID in aaIDs
		]
		)



if __name__ == '__main__':
	main()
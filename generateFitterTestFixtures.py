'''

Let's replace this file with something better - John

'''

print 'Warning!  This will take forever to run.  Get it running on the cluster!'

import cPickle
import os

import wholecell.util.Fitter as wcFitter
import wholecell.sim.Simulation as wcSimulation
import wholecell.sim.logger.Disk as wcDisk

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
nSeeds = 100

logDir = os.path.join('out', 'tests', 'rnaProduction')

kb = cPickle.load(open(KB_PATH, "rb"))

for seed in xrange(nSeeds):
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

print 'Finished'

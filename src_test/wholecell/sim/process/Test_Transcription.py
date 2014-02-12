"""
Test Transcription.py
Tests Transcription process

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/16/2013
"""

import BaseLargeTest
import os

import nose.plugins.attrib as noseAttrib
import numpy
import tables
import scipy.stats

from wholecell.util.Constants import Constants

# TODO: check rRNA production
# TODO: check increase in RNA mass instead of counts

class Test_Transcription(BaseLargeTest.BaseLargeTest):

	# Fixture options for BaseLargeTest class to use
	initRnapCnts = 1000.
	initRnaseCnts = 1000.
	initEnzCnts = 2000.
	ntpCounts = 1e6
	h2oCounts = 1e6
	fixture_opts = {'fixture_dir'		: os.path.join('out', 'test', 'Test_Transcription'),
					'processes' 		: ["Transcription"],
					'free_molecules'	: [["ATP[c]", ntpCounts],
											["UTP[c]", ntpCounts],
											["CTP[c]", ntpCounts],
											["GTP[c]", ntpCounts],
											["EG10893-MONOMER[c]", initEnzCnts],
											["RPOB-MONOMER[c]", initEnzCnts],
											["RPOC-MONOMER[c]", initEnzCnts],
											["RPOD-MONOMER[c]", initEnzCnts]],
					'length_sec'		: 21,
					'n_seeds'			: 1}

	FIXTURE_DIR = fixture_opts['fixture_dir']

	@classmethod
	def setUpClass(cls):

		super(Test_Transcription, cls).setUpClass()

		doublingTime = 3600.
		simLength = 500

		dirs = os.walk(FIXTURE_DIR).next()[1] # get all of the directories

		nSims = len(dirs)

		ntpIDs = ['ATP', 'UTP', 'CTP', 'GTP']
		rrnaIDs = ['RRLA-RRNA:nascent', 'RRLB-RRNA:nascent',
			'RRLC-RRNA:nascent', 'RRLD-RRNA:nascent', 'RRLE-RRNA:nascent',
			'RRLG-RRNA:nascent', 'RRLH-RRNA:nascent']

		assignedIdxs = False

		ntpIdxs = None
		rnaIdxs = None
		rrnaIdxs = None
		processIdx = None
		compartmentIdx = None

		width = 10

		cls.ntpInitialUsage = numpy.zeros((len(ntpIDs), nSims), float)
		cls.ntpFinalUsage = numpy.zeros((len(ntpIDs), nSims), float)

		cls.rnaInitialProduction = None
		cls.rnaFinalProduction = None

		cls.expectedRatio = numpy.exp(numpy.log(2)/doublingTime * simLength)

		for iSim, simDir in enumerate(dirs):
			with tables.openFile(os.path.join(FIXTURE_DIR, simDir, 'MoleculeCounts.hdf')) as h5file:
				if not assignedIdxs:
					names = h5file.get_node('/names')
					indexes = h5file.get_node('/indexes')

					molIDs = names.molIDs.read()
					processes = names.processes.read()
					compartments = names.compartments.read()
					
					ntpIdxs = numpy.array([molIDs.index(id_) for id_ in ntpIDs])
					rnaIdxs = indexes.nascentRnas.read()
					rrnaIdxs = numpy.array([molIDs.index(id_) for id_ in rrnaIDs])

					processIdx = processes.index('Transcription')
					compartmentIdx = compartments.index('c')

					cls.rnaInitialProduction = numpy.zeros((len(rnaIdxs), nSims), float)
					cls.rnaFinalProduction = numpy.zeros((len(rnaIdxs), nSims), float)

					cls.rnaFinalCount = numpy.zeros((len(rnaIdxs), nSims), float)

					assignedIdxs = True

				mc = h5file.root.MoleculeCounts

				cls.ntpInitialUsage[:, iSim] = (mc.read(1, 1+width, None, 'countsBulkPartitioned')[:, ntpIdxs, compartmentIdx, processIdx].mean(0)
					- mc.read(1, width+1, None, 'countsBulkReturned')[:, ntpIdxs, compartmentIdx, processIdx].mean(0))
				cls.ntpFinalUsage[:, iSim] = (mc.read(simLength+1 - width, simLength+1, None, 'countsBulkPartitioned')[:, ntpIdxs, compartmentIdx, processIdx].mean(0)
					- mc.read(simLength+1 - width, simLength+1, None, 'countsBulkReturned')[:, ntpIdxs, compartmentIdx, processIdx].mean(0))

				cls.rnaInitialProduction[:, iSim] = (mc.read(1, 1+width, None, 'countsBulkReturned')[:, rnaIdxs, compartmentIdx, processIdx].mean(0)
					- mc.read(1, width+1, None, 'countsBulkPartitioned')[:, rnaIdxs, compartmentIdx, processIdx].mean(0))
				cls.rnaFinalProduction[:, iSim] = (mc.read(simLength+1 - width, simLength+1, None, 'countsBulkReturned')[:, rnaIdxs, compartmentIdx, processIdx].mean(0)
					- mc.read(simLength+1 - width, simLength+1, None, 'countsBulkPartitioned')[:, rnaIdxs, compartmentIdx, processIdx].mean(0))

				cls.rnaFinalCount[:, iSim] = mc[-1]['countsBulk'][rnaIdxs, compartmentIdx]

				# indexing order for mc.read(...): time, molecule, compartment, partition

		assignedFittedParameters = False

		cls.rnaSynthProb = None

		for iSim, simDir in enumerate(dirs):
			with tables.openFile(os.path.join(FIXTURE_DIR, simDir, 'Main.hdf')) as h5file:
				if not assignedFittedParameters:
					cls.rnaSynthProb = h5file.get_node('/fitParameters').rnaSynthProb.read()

					assignedFittedParameters = True


	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		pass


	def tearDown(self):
		pass


	@noseAttrib.attr('largetest', 'modelfitting', 'transcription')
	def test_ntpUsage(self):
		# Test for exponential NTP usage
		# NOTE: a proper fit should actually exceed the expected ratio since this doesn't include degradation
		ntpRatio = self.ntpFinalUsage/self.ntpInitialUsage
		self.assertTrue(numpy.allclose(self.expectedRatio, ntpRatio, rtol = 0.25))
		self.assertTrue((self.expectedRatio <= ntpRatio).all(), 'Final NTP usage was less than expected.')


	@noseAttrib.attr('largetest', 'modelfitting', 'transcription')
	def test_rnaProduction(self):
		# Test for exponential RNA production
		# NOTE: a proper fit should actually exceed the expected ratio since this doesn't include degradation
		rnaRatio = self.rnaFinalProduction.sum()/self.rnaInitialProduction.sum()
		self.assertTrue(numpy.allclose(self.expectedRatio, rnaRatio, rtol = 0.25))
		self.assertTrue((self.expectedRatio <= rnaRatio).all(), 'Final RNA production was less than expected.')

		# Test for appropriate ratios of rRNA production
		averageCounts = self.rnaFinalCount.mean(1)

		expectedProduction = self.rnaSynthProb * averageCounts.sum()

		# Assert at least a 95% positive correlation
		self.assertGreater(
			scipy.stats.pearsonr(averageCounts, expectedProduction)[0],
			0.95
			)

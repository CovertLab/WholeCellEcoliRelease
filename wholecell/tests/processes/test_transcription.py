"""
Test Transcription.py
Tests Transcription process

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 8/16/2013
"""

from __future__ import division

import wholecell.tests.processes.base_large_test as base_large_test # TODO: better location for this file
import os

import nose.plugins.attrib as noseAttrib
import numpy as np
import tables
import scipy.stats

# TODO: check rRNA production
# TODO: check increase in RNA mass instead of counts

class Test_Transcription(base_large_test.BaseLargeTest):

	# Fixture options for BaseLargeTest class to use
	initRnapCnts = 1000.
	initRnaseCnts = 1000.
	initEnzCnts = 2000.
	ntpCounts = 1e6
	h2oCounts = 1e6
	fixtureOpts = {
		'outputDir':os.path.join('out', 'test', 'Test_Transcription'),
		'includedProcesses':["Transcription"],
		'freeMolecules':[
			["ATP[c]", ntpCounts],
			["UTP[c]", ntpCounts],
			["CTP[c]", ntpCounts],
			["GTP[c]", ntpCounts],
			["EG10893-MONOMER[c]", initEnzCnts],
			["RPOB-MONOMER[c]", initEnzCnts],
			["RPOC-MONOMER[c]", initEnzCnts],
			["RPOD-MONOMER[c]", initEnzCnts]
			],
		'lengthSec':500,
		'nSeeds':8
		}

	fixtureDir = fixtureOpts['outputDir']

	@classmethod
	def setUpClass(cls):
		super(Test_Transcription, cls).setUpClass()

		doublingTime = 3600.
		simLength = cls.fixtureOpts['lengthSec']

		dirs = os.walk(Test_Transcription.fixtureDir).next()[1] # get all of the directories
		nSims = len(dirs)

		# TODO: get these two lines working
		# nSims = cls.fixtureOpts['nSeeds']
		# dirs = [os.path.join(cls.fixtureDir, 'sim{}'.format(i)) for i in range(nSims)]

		ntpIDs = ['ATP', 'UTP', 'CTP', 'GTP']
		rrnaIDs = ['RRLA-RRNA', 'RRLB-RRNA',
			'RRLC-RRNA', 'RRLD-RRNA', 'RRLE-RRNA',
			'RRLG-RRNA', 'RRLH-RRNA']

		assignedIdxs = False

		ntpIdxs = None
		rnaIdxs = None
		rrnaIdxs = None
		processIdx = None
		compartmentIdx = None

		startAt = 2
		width = 10

		cls.ntpInitialUsage = np.zeros((len(ntpIDs), nSims), float)
		cls.ntpFinalUsage = np.zeros((len(ntpIDs), nSims), float)

		cls.rnaInitialProduction = None
		cls.rnaFinalProduction = None

		cls.expectedRatio = np.exp(np.log(2)/doublingTime * simLength)

		for iSim, simDir in enumerate(dirs):
			with tables.openFile(os.path.join(cls.fixtureDir, simDir, 'MoleculeCounts.hdf')) as h5file:
				if not assignedIdxs:
					names = h5file.get_node('/names')
					indexes = h5file.get_node('/indexes')

					molIDs = names.molIDs.read()
					processes = names.processes.read()
					compartments = names.compartments.read()
					
					ntpIdxs = np.array([molIDs.index(id_) for id_ in ntpIDs])
					rnaIdxs = indexes.rnas.read()
					rrnaIdxs = np.array([molIDs.index(id_) for id_ in rrnaIDs])

					processIdx = processes.index('Transcription')
					compartmentIdx = compartments.index('c')

					cls.rnaInitialProduction = np.zeros((len(rnaIdxs), nSims), float)
					cls.rnaFinalProduction = np.zeros((len(rnaIdxs), nSims), float)

					cls.rnaFinalCount = np.zeros((len(rnaIdxs), nSims), float)

					assignedIdxs = True

				mc = h5file.root.MoleculeCounts

				cls.ntpInitialUsage[:, iSim] = (mc.read(startAt, startAt+width, None, 'countsBulkPartitioned')[:, ntpIdxs, compartmentIdx, processIdx].mean(0)
					- mc.read(startAt, startAt+width, None, 'countsBulkReturned')[:, ntpIdxs, compartmentIdx, processIdx].mean(0))
				cls.ntpFinalUsage[:, iSim] = (mc.read(simLength+1 - width, simLength+1, None, 'countsBulkPartitioned')[:, ntpIdxs, compartmentIdx, processIdx].mean(0)
					- mc.read(simLength+1 - width, simLength+1, None, 'countsBulkReturned')[:, ntpIdxs, compartmentIdx, processIdx].mean(0))

				cls.rnaInitialProduction[:, iSim] = (mc.read(startAt, startAt+width, None, 'countsBulkReturned')[:, rnaIdxs, compartmentIdx, processIdx].mean(0)
					- mc.read(startAt, startAt+width, None, 'countsBulkPartitioned')[:, rnaIdxs, compartmentIdx, processIdx].mean(0))
				cls.rnaFinalProduction[:, iSim] = (mc.read(simLength+1 - width, simLength+1, None, 'countsBulkReturned')[:, rnaIdxs, compartmentIdx, processIdx].mean(0)
					- mc.read(simLength+1 - width, simLength+1, None, 'countsBulkPartitioned')[:, rnaIdxs, compartmentIdx, processIdx].mean(0))

				cls.rnaFinalCount[:, iSim] = mc[-1]['countsBulk'][rnaIdxs, compartmentIdx]

				# indexing order for mc.read(...): time, molecule, compartment, partition

		assignedFittedParameters = False

		cls.rnaSynthProb = None

		for iSim, simDir in enumerate(dirs):
			with tables.openFile(os.path.join(cls.fixtureDir, simDir, 'Main.hdf')) as h5file:
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
		ntpRatio = self.ntpFinalUsage/self.ntpInitialUsage

		self.assertTrue(np.allclose(self.expectedRatio, ntpRatio, rtol = 0.05))


	@noseAttrib.attr('largetest', 'modelfitting', 'transcription')
	def test_rnaProduction(self):
		# Test for exponential RNA production
		rnaRatio = self.rnaFinalProduction.sum()/self.rnaInitialProduction.sum()
		self.assertTrue(np.allclose(self.expectedRatio, rnaRatio, rtol = 0.05))

		# Test for appropriate ratios of rRNA production
		averageCounts = self.rnaFinalCount.mean(1)

		expectedProduction = self.rnaSynthProb * averageCounts.sum()

		# Assert at least a 95% positive correlation
		self.assertGreater(
			scipy.stats.pearsonr(averageCounts, expectedProduction)[0],
			0.95
			)

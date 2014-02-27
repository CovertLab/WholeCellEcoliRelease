#!/usr/bin/env python

"""
Test reduced simulations.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/6/20
"""

import os
import cPickle
import unittest

import nose.plugins.attrib as noseAttrib

import wholecell.utils.fitter as wcFitter
import wholecell.sim.simulation as wcSimulation

KB_PATH = os.path.join('data', 'fixtures', 'KnowledgeBase.cPickle')

class Test_reducedSimulations(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		pass

	def tearDown(self):
		pass

	@noseAttrib.attr('mediumtest')
	def test_transcriptionSimulation(self):
		sim = wcSimulation.Simulation(
			includedProcesses = ["Transcription"],
			freeMolecules = [
				["ATP[c]", 1e6],
				["UTP[c]", 1e6],
				["CTP[c]", 1e6],
				["GTP[c]", 1e6],
				["EG10893-MONOMER[c]", 2000.0],
				["RPOB-MONOMER[c]", 2000.0],
				["RPOC-MONOMER[c]", 2000.0],
				["RPOD-MONOMER[c]", 2000.0]
			],
			lengthSec = 10
			)
		
		self.assertEqual(
			sim.processes.viewkeys(),
			set(['Transcription', 'FreeProduction'])
			)

		sim.run()

		mc = sim.states['MoleculeCounts']

		ntpView = mc.countsBulkViewNew(
			["ATP[c]", "UTP[c]", "CTP[c]", "GTP[c]"]
			)

		ntpMin = 0.5e6 # allow some tolerance for process usage

		self.assertTrue(
			(ntpView.countsBulk() > ntpMin).all()
			)

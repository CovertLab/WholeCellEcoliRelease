#!/usr/bin/env python

"""
Test reduced simulations.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/6/20
"""

from __future__ import division

import os
import cPickle
import unittest

import nose.plugins.attrib as noseAttrib

import wholecell.sim.simulation as wcSimulation

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

	@noseAttrib.attr('mediumtest', 'reducedSimulation')
	def test_transcriptionSimulation(self):
		ntpLevels = 1e6
		enzLevels = 2000.

		sim = wcSimulation.Simulation(
			includedProcesses = ["Transcription"],
			freeMolecules = [
				["ATP[c]", ntpLevels],
				["UTP[c]", ntpLevels],
				["CTP[c]", ntpLevels],
				["GTP[c]", ntpLevels],
				["EG10893-MONOMER[c]", enzLevels],
				["RPOB-MONOMER[c]", enzLevels],
				["RPOC-MONOMER[c]", enzLevels],
				["RPOD-MONOMER[c]", enzLevels]
			],
			lengthSec = 10,
			reconstructKB = True,
			seed = 1
			)
		
		self.assertEqual(
			sim.processes.viewkeys(),
			set(['Transcription', 'FreeProduction'])
			)

		sim.run()

		bulkMolecules = sim.states['BulkMolecules']

		ntpView = sim.processes['Transcription'].ntps

		ntpMin = ntpLevels/2 # allow some tolerance for process usage

		self.assertTrue(
			(ntpView.counts() > ntpMin).all()
			)

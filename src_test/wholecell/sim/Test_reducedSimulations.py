#!/usr/bin/env python

"""
Test reduced simulations.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/6/20
"""

import unittest

import wholecell.util.SimulationParser as wcSimulationParser

import nose.plugins.attrib as noseAttrib

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

	@noseAttrib.attr('working')
	@noseAttrib.attr('smalltest')
	def test_transcriptionSimulation(self):
		sim = wcSimulationParser.parseSimulationFromJsonString(
			'''
			{
				"processesToInclude":["Transcription"],
				"useShellLogger":false,
				"fitSimulation":true,
				"simOptions":{
					"seed":10,
					"lengthSec":10
				},
				"freeMolecules":[
					["ATP[c]", 1e6],
					["UTP[c]", 1e6],
					["CTP[c]", 1e6],
					["GTP[c]", 1e6],
					["EG10893-MONOMER[c]", 2000.0],
					["RPOB-MONOMER[c]", 2000.0],
					["RPOC-MONOMER[c]", 2000.0],
					["RPOD-MONOMER[c]", 2000.0]
				]
			}
			'''
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

#!/usr/bin/env python

"""
Test SimulationParser.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/6/2014
"""

import unittest

import wholecell.sim.simulation as wcSimulation
import wholecell.utils.simulation_parser as wcSimulationParser

import nose.plugins.attrib as noseAttrib

class Test_SimulationParser(unittest.TestCase):

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

	
	@noseAttrib.attr('smalltest')
	def test_defaults(self):
		sim = wcSimulationParser.parseSimulationFromJsonString('{}')

		# TODO: assert fitting
		# TODO: assert useCachedKB/cachedKB
		# TODO: assert shell logger
		# TODO: assert not disk logger
		# TODO: assert no sim options

		self.assertEqual(
			sim.processes.viewkeys(),
			set(wcSimulation.DEFAULT_PROCESSES)
			)

		self.assertEqual(
			sim.freeMolecules,
			None
			)

	
	@noseAttrib.attr('smalltest')
	def test_reducedSimulation(self):
		sim = wcSimulationParser.parseSimulationFromJsonString(
			'''
			{
				"processesToInclude":["Transcription"],
				"useShellLogger":false,
				"fitSimulation":false,
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

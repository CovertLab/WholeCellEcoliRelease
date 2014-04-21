"""
Test Simulation.py
Tests whole-cell simulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/8/2013
"""

from __future__ import division

import unittest
import warnings
import nose.plugins.attrib as noseAttrib

import numpy as np
import cPickle
import os

import wholecell.sim.simulation
import wholecell.loggers.disk
import wholecell.loggers.shell

import wholecell.utils.config
TEST_FIXTURE_DIR = wholecell.utils.config.TEST_FIXTURE_DIR

# TODO: add "short sim" fixture

class Test_Simulation(unittest.TestCase):
	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.sim = cPickle.load(open(os.path.join(TEST_FIXTURE_DIR, "Simulation.cPickle"), "r"))
		import wholecell.utils.knowledgebase_fixture_manager
		self.kb = wholecell.utils.knowledgebase_fixture_manager.loadKnowledgeBase(os.path.join(TEST_FIXTURE_DIR, 'KnowledgeBase.cPickle'))
		
	def tearDown(self):
		pass


	# --- Tests for run-time errors ---
	@noseAttrib.attr('mediumtest', 'simulation')
	def test_construction(self):
		# Construct simulation
		sim = wholecell.sim.simulation.Simulation(reconstructKB = True)


	@noseAttrib.attr('mediumtest', 'simulation')
	def test_run(self):
		# Simulate
		sim = wholecell.sim.simulation.Simulation(seed = 0, lengthSec = 10, reconstructKB = True)
		sim.run()

		self.assertEqual(10, sim.time())


	@noseAttrib.attr('mediumtest', 'saveload', 'simulation')
	def test_disk_logger(self):
		# Output directory
		outDir = os.path.join("out", "test", "SimulationTest_testLogging")

		# Run simulation
		sim = wholecell.sim.simulation.Simulation(
			seed = 0, lengthSec = 10, logToDisk = True, outputDir = outDir,
			overwriteExistingFiles = True,
			reconstructKB = True,
			includedStates = ['Mass', 'BulkMolecules', 'UniqueMolecules', 'Chromosome', 'Transcripts']
			)

		sim.run()
		
		reloadedSim = wholecell.sim.simulation.Simulation.loadSimulation(outDir, timePoint = 10)

		self.assertEqual(
			sim.states['BulkMolecules'].container,
			reloadedSim.states['BulkMolecules'].container,
			)

		self.assertEqual(
			sim.states['UniqueMolecules'].container,
			reloadedSim.states['UniqueMolecules'].container,
			)

		self.assertEqual(
			sim.states['Transcripts'].container,
			reloadedSim.states['Transcripts'].container,
			)

		self.assertEqual(
			sim.states['Chromosome'].container,
			reloadedSim.states['Chromosome'].container,
			)


	# --- Test ability to remove processes from simulation ---
	@noseAttrib.attr('mediumtest')
	def test_removeProcesses(self):
		sim = wholecell.sim.simulation.Simulation(includedProcesses = ['Transcription'], reconstructKB = True)
		self.assertEqual(['Transcription'], sim.processes.keys())

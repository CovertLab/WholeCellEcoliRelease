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


	@noseAttrib.attr('mediumtest', 'simulation')
	def test_removeProcesses(self):
		# Test ability to remove processes from simulation
		sim = wholecell.sim.simulation.Simulation(includedProcesses = ['BulkTranscription'], reconstructKB = True)
		self.assertEqual(['BulkTranscription'], sim.processes.keys())

	
	@noseAttrib.attr('mediumtest', 'simulation')
	def test_removeStates(self):
		# Test ability to remove states from simulation
		sim = wholecell.sim.simulation.Simulation(includedProcesses = ['BulkTranscription'], includedStates = ['BulkMolecules'], reconstructKB = True)
		self.assertEqual(['BulkMolecules'], sim.states.keys())

	
	@noseAttrib.attr('mediumtest', 'simulation')
	def test_processInclusionOrder(self):
		# Test that included processes follow the proscribed order (important for other tests)
		sim = wholecell.sim.simulation.Simulation(includedProcesses = ['BulkTranscription', 'Translation'], reconstructKB = True)
		self.assertEqual(['BulkTranscription', 'Translation'], sim.processes.keys())

		sim = wholecell.sim.simulation.Simulation(includedProcesses = ['Translation', 'BulkTranscription'], reconstructKB = True)
		self.assertEqual(['Translation', 'BulkTranscription'], sim.processes.keys())


	
	@noseAttrib.attr('mediumtest', 'simulation')
	def test_processEvaluationOrderIndependence(self):
		# Test that process evaluation order does not determine the consequences of one time step
		# NOTE: comparing independent states beyond BulkMolecules is a tricky affair

		from wholecell.sim.simulation import DEFAULT_PROCESSES

		processes1 = DEFAULT_PROCESSES[:]
		processes2 = DEFAULT_PROCESSES[::-1]

		sim1 = wholecell.sim.simulation.Simulation(
			includedProcesses = processes1,
			reconstructKB = True,
			seed = 0,
			lengthSec = 10
			)

		sim1.run()

		sim1Mass = sim1.states['Mass']

		masses1 = [
			sim1Mass.metabolite,
			sim1Mass.rna,
			sim1Mass.protein
			]

		del sim1

		sim2 = wholecell.sim.simulation.Simulation(
			includedProcesses = processes2,
			reconstructKB = True,
			seed = 0,
			lengthSec = 10
			)

		sim2.run()

		sim2Mass = sim2.states['Mass']

		masses2 = [
			sim2Mass.metabolite,
			sim2Mass.rna,
			sim2Mass.protein
			]

		del sim2

		self.assertEqual(
			masses1,
			masses2
			)






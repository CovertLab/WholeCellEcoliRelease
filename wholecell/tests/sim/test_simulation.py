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

# TODO: add "short sim" fixture

class Test_Simulation(unittest.TestCase):
	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.sim = cPickle.load(open(os.path.join("fixtures", "Simulation.cPickle"), "r"))
		self.kb = cPickle.load(open(os.path.join("fixtures","KnowledgeBase.cPickle"), "r"))

	def tearDown(self):
		pass


	# --- Tests for run-time errors ---
	@noseAttrib.attr('smalltest')
	def test_construction(self):
		# Construct simulation
		sim = wholecell.sim.simulation.Simulation()

	@noseAttrib.attr('mediumtest')
	def test_run(self):
		# Simulate
		sim = wholecell.sim.simulation.Simulation(seed = 0, lengthSec = 10)
		sim.run()

		self.assertEqual(10, sim.states["Time"].value)

	@noseAttrib.attr('mediumtest')
	def test_disk_logger(self): #_and_shell_logger(self):
		# Output directory
		outDir = os.path.join("out", "test", "SimulationTest_testLogging")

		# Run simulation
		sim = wholecell.sim.simulation.Simulation(
			seed = 0, lengthSec = 10, logToDisk = True, outputDir = outDir,
			overwriteExistingFiles = True
			)
		sim.run()
		
		reloadedSim = wholecell.sim.simulation.Simulation.loadSimulation(outDir, timePoint = 10)

		state_keys = sim.states.keys()
		# Need to check RandStream in another way
		state_keys.pop(state_keys.index('RandStream'))
		for state_id in state_keys:
			dynamics_keys = sim.states[state_id].getDynamics().keys()

			if 'growth' in dynamics_keys:
				# Growth calculated based on difference in time-steps will not match up
				dynamics_keys.pop(dynamics_keys.index('growth'))

			for d_key in dynamics_keys:
				if isinstance(sim.states[state_id].getDynamics()[d_key], np.ndarray):
					self.assertEqual(sim.states[state_id].getDynamics()[d_key].tolist(),
						reloadedSim.states[state_id].getDynamics()[d_key].tolist())
				else:
					self.assertEqual(sim.states[state_id].getDynamics()[d_key],
						reloadedSim.states[state_id].getDynamics()[d_key])
		# Check RandStream
		self.assertEqual(sim.states['RandStream'].getDynamics()['value'][1].tolist(),
						reloadedSim.states['RandStream'].getDynamics()['value'][1].tolist())


	@noseAttrib.attr('mediumtest')
	def test_reload_at_later_timepoint(self):
		# Output directory
		outDir = os.path.join("out", "test", "SimulationTest_test_reload_at_later_timepoint")

		lengthSec = 10.

		# Run simulation
		sim = wholecell.sim.simulation.Simulation(
			seed = 0, lengthSec = lengthSec, logToDisk = True, outputDir = outDir,
			overwriteExistingFiles = True
			)
		sim.run()

		# TODO: Finish - call from Simulation.Simulation.loadSimulation
		reloadedSim = wholecell.sim.simulation.Simulation.loadSimulation(outDir, timePoint = 5)

		self.assertEqual(reloadedSim.lengthSec, lengthSec)

		self.assertEqual(reloadedSim.initialStep, 5)
		self.assertEqual(reloadedSim.states['Time'].value, 5.)
		reloadedSim.run()

		state_keys = sim.states.keys()
		# Need to check RandStream in another way
		state_keys.pop(state_keys.index('RandStream'))
		for state_id in state_keys:
			dynamics_keys = sim.states[state_id].getDynamics().keys()

			if 'growth' in dynamics_keys:
				# Growth calculated based on difference in time-steps will not match up
				dynamics_keys.pop(dynamics_keys.index('growth'))

			for d_key in dynamics_keys:
				if isinstance(sim.states[state_id].getDynamics()[d_key], np.ndarray):
					self.assertEqual(sim.states[state_id].getDynamics()[d_key].tolist(),
						reloadedSim.states[state_id].getDynamics()[d_key].tolist())
				else:
					self.assertEqual(sim.states[state_id].getDynamics()[d_key],
						reloadedSim.states[state_id].getDynamics()[d_key])
		# Check RandStream
		self.assertEqual(sim.states['RandStream'].getDynamics()['value'][1].tolist(),
						reloadedSim.states['RandStream'].getDynamics()['value'][1].tolist())


	@noseAttrib.attr('smalltest')
	def test_getDynamics(self):
		sim = self.sim
		dynamics = sim.getDynamics()
		self.assertEqual(
			dynamics.viewkeys(),
			{'RandStream', 'UniqueMolecules', 'Mass', 'BulkMolecules', 'Time', 'Chromosome'}
			)


	# --- Test ability to remove processes from simulation ---
	@noseAttrib.attr('smalltest')
	def test_removeProcesses(self):
		sim = wholecell.sim.simulation.Simulation(includedProcesses = ['Transcription'])
		self.assertEqual(['Transcription'], sim.processes.keys())

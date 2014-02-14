"""
Test Simulation.py
Tests whole-cell simulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/8/2013
"""

import unittest
import warnings
import nose.plugins.attrib as noseAttrib

import numpy
import cPickle
import os

import wholecell.sim.simulation
import wholecell.loggers.disk
import wholecell.loggers.shell

class Test_Simulation(unittest.TestCase):
	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		self.sim = cPickle.load(open(os.path.join("data", "fixtures", "Simulation.cPickle"), "r"))
		self.kb = cPickle.load(open(os.path.join("data","fixtures","KnowledgeBase.cPickle"), "r"))

	def tearDown(self):
		pass


	# --- Tests for run-time errors ---
	@noseAttrib.attr('smalltest')
	def test_construction(self):
		# Construct simulation
		sim = wholecell.sim.simulation.Simulation()

	@noseAttrib.attr('smalltest')
	def test_initialize(self):
		# Construct simulation
		sim = wholecell.sim.simulation.Simulation()
		sim.initialize(self.kb)

	@noseAttrib.attr('mediumtest')
	def test_run(self):
		# Simulate
		sim = self.sim
		sim.initialize(self.kb)
		sim.setOptions({"lengthSec": 10})
		sim.run()

		self.assertEqual(10, sim.states["Time"].value)

	@noseAttrib.attr('mediumtest')
	def test_disk_logger(self): #_and_shell_logger(self):
		# Output directory
		outDir = os.path.join("out", "test", "SimulationTest_testLogging")

		# Run simulation
		sim = self.sim
		sim.initialize(self.kb)
		sim.setOptions({"lengthSec": 10})
		sim.loggerAdd(wholecell.loggers.disk.Disk(outDir = outDir, allowOverwrite = True))
		sim.run()
		
		# TODO: Finish - call from Simulation.Simulation.loadSimulation
		reloadedSim = wholecell.sim.simulation.Simulation.loadSimulation(self.kb, outDir, timePoint = 10)

		state_keys = sim.states.keys()
		# Need to check RandStream in another way
		state_keys.pop(state_keys.index('RandStream'))
		for state_id in state_keys:
			dynamics_keys = sim.states[state_id].getDynamics().keys()

			if 'growth' in dynamics_keys:
				# Growth calculated based on difference in time-steps will not match up
				dynamics_keys.pop(dynamics_keys.index('growth'))

			for d_key in dynamics_keys:
				if isinstance(sim.states[state_id].getDynamics()[d_key], numpy.ndarray):
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

		# Run simulation
		sim = self.sim
		sim.initialize(self.kb)
		sim.setOptions({"lengthSec": 10})
		sim.loggerAdd(wholecell.loggers.shell.Shell())
		sim.loggerAdd(wholecell.loggers.disk.Disk(outDir = outDir, allowOverwrite = True))
		sim.run()

		# TODO: Finish - call from Simulation.Simulation.loadSimulation
		reloadedSim = wholecell.sim.simulation.Simulation.loadSimulation(self.kb, outDir, timePoint = 5)
		reloadedSim.setOptions({"lengthSec": 10})
		reloadedSim.loggerAdd(wholecell.loggers.shell.Shell())

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
				if isinstance(sim.states[state_id].getDynamics()[d_key], numpy.ndarray):
					self.assertEqual(sim.states[state_id].getDynamics()[d_key].tolist(),
						reloadedSim.states[state_id].getDynamics()[d_key].tolist())
				else:
					self.assertEqual(sim.states[state_id].getDynamics()[d_key],
						reloadedSim.states[state_id].getDynamics()[d_key])
		# Check RandStream
		self.assertEqual(sim.states['RandStream'].getDynamics()['value'][1].tolist(),
						reloadedSim.states['RandStream'].getDynamics()['value'][1].tolist())

	@noseAttrib.attr('mediumtest')
	def test_loadSimulation_method(self):
		with self.assertRaises(Exception) as context:
			sim = self.sim
			sim.initialize(self.kb)
			sim.setOptions({"lengthSec": 2})
			outDir = os.path.join("out", "test", "SimulationTest_testLogging")
			sim.loggerAdd(wholecell.loggers.disk.Disk(outDir = outDir, allowOverwrite = True))
			sim.run()
			wholecell.sim.simulation.Simulation.loadSimulation(self.kb, outDir, timePoint = 3)

		self.assertEqual(context.exception.message, 'Time point chosen to load is out of range!\n')

	@noseAttrib.attr('smalltest')
	def test_get_and_set_options(self):
		sim = self.sim
		options = sim.getOptions()
		self.assertEqual(options.keys(), ['states', 'processes', 'seed', 'lengthSec', 'timeStepSec'])
		sim.setOptions({'lengthSec' : 9., 'timeStepSec' : 3.})
		options = sim.getOptions()
		self.assertEqual(options['lengthSec'], 9.)
		self.assertEqual(options['timeStepSec'], 3.)

		with self.assertRaises(Exception) as context:
			sim.setOptions({'fubar' : True})
		self.assertEqual(context.exception.message, "Invalid options:\n -%s" % "fubar")

		# Just making sure there is no error thrown in the Simulation.py code.
		sim.setOptions({'processes' : {'Metabolism' : {'lpSolver' : 'test'}}})
		self.assertEqual(sim.processes['Metabolism'].lpSolver, 'test')

		sim.states['Time'].meta['options'] = ['test']
		sim.setOptions({'states' : {'Time' : {'test' : 'test_val'}}})
		self.assertEqual(sim.states['Time'].test, 'test_val')

	@noseAttrib.attr('smalltest')
	def test_getDynamics(self):
		sim = self.sim
		dynamics = sim.getDynamics()
		self.assertEqual(dynamics.keys(), ['RandStream', 'Mass', 'MoleculeCounts', 'Time'])


	# --- Test ability to remove processes from simulation ---
	@noseAttrib.attr('smalltest')
	def test_removeProcesses(self):
		sim = wholecell.sim.simulation.Simulation(processesToInclude = ['Transcription'])
		sim.initialize(self.kb)
		self.assertEqual(['Transcription'], sim.processes.keys())

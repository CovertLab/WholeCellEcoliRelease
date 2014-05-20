#!/usr/bin/env python

"""
Simulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

from __future__ import division

import collections
import os

import numpy as np
import tables

import wholecell.utils.rand_stream
import wholecell.utils.config
import wholecell.utils.constants
import wholecell.utils.knowledgebase_fixture_manager
import wholecell.reconstruction.fitter
import wholecell.sim.sim_definition


class SimulationException(Exception):
	pass


class Simulation(object):
	""" Simulation """

	# Constructors
	def __init__(self, simDefinition = None, **kwargs):
		# Establish simulation options
		if simDefinition is not None and kwargs:
			raise SimulationException(
				"Simulations cannot be instantiated with both a SimDefinition instance and keyword arguments"
				)

		elif simDefinition is None:
			simDefinition = wholecell.sim.sim_definition.SimDefinition(**kwargs)

		self._options = simDefinition

		# Set time parameters
		self.lengthSec = self._options.lengthSec
		self.timeStepSec = self._options.timeStepSec

		self.initialStep = 0
		self.simulationStep = 0

		# Set random seed
		if self._options.seed is None:
			# This is roughly consistent with how NumPy choose a seed if none
			# is provided
			import random
			self.seed = random.SystemRandom().randrange(2**32-1)

		else:
			self.seed = self._options.seed

		self.randStream = wholecell.utils.rand_stream.RandStream(seed = self.seed)

		# Create KB
		self.kbDir = wholecell.utils.constants.SIM_FIXTURE_DIR

		if not self._options.rebuildKB:
			kb = wholecell.utils.knowledgebase_fixture_manager.loadKnowledgeBase(os.path.join(self.kbDir, 'KnowledgeBase.cPickle'))
		else:
			kb = wholecell.utils.knowledgebase_fixture_manager.cacheKnowledgeBase(self.kbDir)

		# Fit KB parameters
		wholecell.reconstruction.fitter.fitKb(kb)

		# Initialize simulation from fit KB
		self._initialize(kb)


	@classmethod
	def initFromFile(cls, filePath, **kwargs):
		import json

		try:
			jsonArgs = json.load(open(filePath))

		except ValueError:
			raise Exception('Caught ValueError; these can be caused by excess commas in the json file, which may not be caught by the syntax checker in your text editor.')

		jsonArgs.update(kwargs)

		return cls(**jsonArgs)


	# Link states and processes
	def _initialize(self, kb):
		self.states = self._options.createStates()
		self.processes = self._options.createProcesses()
		self.listeners = self._options.createListeners()
		self.hooks = self._options.createHooks()

		for state in self.states.itervalues():
			state.initialize(self, kb)

		for process in self.processes.itervalues():
			process.initialize(self, kb)

		for listener in self.listeners.itervalues():
			listener.initialize(self, kb)

		for hook in self.hooks.itervalues():
			hook.initialize(self, kb)

		for state in self.states.itervalues():
			state.allocate()

		for listener in self.listeners.itervalues():
			listener.allocate()
		
		from wholecell.reconstruction.initial_conditions import calcInitialConditions

		calcInitialConditions(self, kb)

		for hook in self.hooks.itervalues():
			hook.postCalcInitialConditions(self)

		# Create loggers
		self.loggers = self._options.createLoggers()

	# -- Run simulation --

	# Run simulation
	def run(self):
		# Perform initial listener update
		for listener in self.listeners.itervalues():
			listener.initialUpdate()

		# Start logging
		for logger in self.loggers:
			logger.initialize(self)

		# Simulate
		while self.time() < self.lengthSec:
			self.simulationStep += 1

			self._evolveState()

		# Run post-simulation hooks
		for hook in self.hooks.itervalues():
			hook.finalize(self)

		# Finish logging
		for logger in self.loggers:
			logger.finalize(self)


	# Calculate temporal evolution
	def _evolveState(self):
		# Update randstreams
		for stateName, state in self.states.iteritems():
			state.seed = self._seedFromName(stateName)
			state.randStream = wholecell.utils.rand_stream.RandStream(
				seed = state.seed
				)

		for processName, process in self.processes.iteritems():
			process.seed = self._seedFromName(processName)
			process.randStream = wholecell.utils.rand_stream.RandStream(
				seed = process.seed
				)

		# TODO: randstreams for hooks?

		# Run pre-evolveState hooks
		for hook in self.hooks.itervalues():
			hook.preEvolveState(self)

		# Update queries
		for state in self.states.itervalues():
			state.updateQueries()

		# Calculate requests
		for process in self.processes.itervalues():
			process.calculateRequest()

		# Partition states among processes
		for state in self.states.itervalues():
			state.partition()

		# Simulate submodels
		for process in self.processes.itervalues():
			process.evolveState()

		# Merge state
		for state in self.states.itervalues():
			state.merge()

		# Update listeners
		for listener in self.listeners.itervalues():
			listener.update()

		# Run post-evolveState hooks
		for hook in self.hooks.itervalues():
			hook.postEvolveState(self)

		# Append loggers
		for logger in self.loggers:
			logger.append(self)


	def _seedFromName(self, name):
		return np.uint32(self.seed + self.simulationStep + hash(name))


	# Save to/load from disk
	def pytablesCreate(self, h5file, expectedRows):
		groupNames = h5file.create_group(
			h5file.root,
			'names',
			'State and process names'
			)

		h5file.create_array(groupNames, 'states', [s for s in self.states.viewkeys()])
		
		if self.processes:
			h5file.create_array(groupNames, 'processes', [s for s in self.processes.viewkeys()])

		# TODO: cache KB


	def pytablesAppend(self, h5file):
		# Included for consistency, eventual features...
		pass


	def pytablesLoad(self, h5file, timePoint):
		pass


	@classmethod
	def loadSimulation(cls, simDir, timePoint, newDir = None, overwriteExistingFiles = False):
		newSim = cls.initFromFile(
			os.path.join(simDir, 'simOpts.json'),
			logToDisk = newDir is not None,
			overwriteExistingFiles = overwriteExistingFiles,
			outputDir = newDir
			)

		with tables.open_file(os.path.join(simDir, 'Main.hdf')) as h5file:
			newSim.pytablesLoad(h5file, timePoint)

		for stateName, state in newSim.states.viewitems():
			with tables.open_file(os.path.join(simDir, stateName + '.hdf')) as h5file:
				state.pytablesLoad(h5file, timePoint)

		for listenerName, listener in newSim.listeners.viewitems():
			with tables.open_file(os.path.join(simDir, listenerName + '.hdf')) as h5file:
				listener.pytablesLoad(h5file, timePoint)

		newSim.initialStep = timePoint

		return newSim


	def options(self):
		return self._options


	def time(self):
		return self.timeStepSec * (self.initialStep + self.simulationStep)


	def timeStep(self):
		return self.initialStep + self.simulationStep

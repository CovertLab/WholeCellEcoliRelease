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
import shutil
import cPickle
import time

import numpy as np
import tables

import wholecell.utils.config
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

		self.randomState = np.random.RandomState(seed = self.seed)

		# Load KB
		kb = cPickle.load(open(self._options.kbLocation, "rb"))

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
		self.loggers = self._options.createLoggers()

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

		# Make permanent reference to evaluation time listener
		self._evalTime = self.listeners["EvaluationTime"]

	# -- Run simulation --

	# Run simulation
	def run(self):
		# Perform initial mass calculations
		for state in self.states.itervalues():
			state.calculatePreEvolveStateMass()
			state.calculatePostEvolveStateMass()

		# Perform initial listener update
		for listener in self.listeners.itervalues():
			listener.initialUpdate()

		# Start logging
		for logger in self.loggers.itervalues():
			logger.initialize(self)

		# Simulate
		while self.time() < self.lengthSec:
			self.simulationStep += 1

			self._evolveState()

		# Run post-simulation hooks
		for hook in self.hooks.itervalues():
			hook.finalize(self)

		# Finish logging
		for logger in self.loggers.itervalues():
			logger.finalize(self)


	# Calculate temporal evolution
	def _evolveState(self):
		# Update randstreams
		for stateName, state in self.states.iteritems():
			state.seed = self._seedFromName(stateName)
			state.randomState = np.random.RandomState(seed = state.seed)

		for processName, process in self.processes.iteritems():
			process.seed = self._seedFromName(processName)
			process.randomState = np.random.RandomState(seed = process.seed)

		# TODO: randstreams for hooks?

		# Run pre-evolveState hooks
		for hook in self.hooks.itervalues():
			hook.preEvolveState(self)

		# Update queries
		# TODO: context manager/function calls for this logic?
		for i, state in enumerate(self.states.itervalues()):
			t = time.time()
			state.updateQueries()
			self._evalTime.updateQueries_times[i] = time.time() - t

		# Calculate requests
		for i, process in enumerate(self.processes.itervalues()):
			t = time.time()
			process.calculateRequest()
			self._evalTime.calculateRequest_times[i] = time.time() - t

		# Partition states among processes
		for i, state in enumerate(self.states.itervalues()):
			t = time.time()
			state.partition()
			self._evalTime.partition_times[i] = time.time() - t

		# Calculate mass of partitioned molecules
		for state in self.states.itervalues():
			state.calculatePreEvolveStateMass()

		# Update listeners
		for listener in self.listeners.itervalues():
			listener.updatePostRequest()

		# Simulate submodels
		for i, process in enumerate(self.processes.itervalues()):
			t = time.time()
			process.evolveState()
			self._evalTime.evolveState_times[i] = time.time() - t

		# Merge state
		for i, state in enumerate(self.states.itervalues()):
			t = time.time()
			state.merge()
			self._evalTime.merge_times[i] = time.time() - t

		# Calculate mass of partitioned molecules, after evolution
		for state in self.states.itervalues():
			state.calculatePostEvolveStateMass()

		# Update listeners
		for listener in self.listeners.itervalues():
			listener.update()

		# Run post-evolveState hooks
		for hook in self.hooks.itervalues():
			hook.postEvolveState(self)

		# Append loggers
		for logger in self.loggers.itervalues():
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

		# Note: the '.encode("ascii")' is necessary because if a json file was parsed
		# for simulation options, it reads things in as unicode (which pytables doesn't like)
		h5file.create_array(groupNames, 'states', [s.encode("ascii") for s in self.states.viewkeys()])
		
		if self.processes:
			h5file.create_array(groupNames, 'processes', [s.encode("ascii") for s in self.processes.viewkeys()])


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

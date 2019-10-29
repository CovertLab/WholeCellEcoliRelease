#!/usr/bin/env python

"""
Simulation

@organization: Covert Lab, Department of Bioengineering, Stanford University
"""

from __future__ import absolute_import
from __future__ import division

import collections
import cPickle
import time

import numpy as np

from wholecell.listeners.evaluation_time import EvaluationTime
from wholecell.utils import filepath

import wholecell.loggers.shell
import wholecell.loggers.disk

DEFAULT_SIMULATION_KWARGS = dict(
	seed = 0,
	lengthSec = 3*60*60, # 3 hours max
	initialTime = 0.,
	massDistribution = True,
	dPeriodDivision = False,
	growthRateNoise = False,
	translationSupply = True,
	variable_elongation_translation = False,
	variable_elongation_transcription = False,
	timeStepSafetyFraction = 1.3,
	maxTimeStep = 0.9,#2.0, # TODO: Reset to 2 once we update PopypeptideElongation
	updateTimeStepFreq = 5,
	logToShell = True,
	logToDisk = False,
	outputDir = None,
	overwriteExistingFiles = False,
	logToDiskEvery = 1,
	simDataLocation = None,
	inheritedStatePath = None,
	)

def _orderedAbstractionReference(iterableOfClasses):
	return collections.OrderedDict(
		(cls.name(), cls())
		for cls in iterableOfClasses
		)


class SimulationException(Exception):
	pass


DEFAULT_LISTENER_CLASSES = (
	EvaluationTime,
	)

class Simulation(object):
	""" Simulation """

	# Attributes that must be set by a subclass
	_definedBySubclass = (
		"_internalStateClasses",
		"_externalStateClasses",
		"_processClasses",
		"_initialConditionsFunction",
		)

	# Attributes that may be optionally overwritten by a subclass
	_listenerClasses = ()
	_hookClasses = ()
	_timeStepSec = .2
	_shellColumnHeaders = ("Time (s)",)

	# Constructors
	def __init__(self, **kwargs):
		# Validate subclassing
		for attrName in self._definedBySubclass:
			if not hasattr(self, attrName):
				raise SimulationException("Simulation subclasses must define"
				+ " the {} attribute.".format(attrName))

		for listenerClass in DEFAULT_LISTENER_CLASSES:
			if listenerClass in self._listenerClasses:
				raise SimulationException("The {} listener is included by"
					+ " default in the Simulation class.".format(
						listenerClass.name())
					)

		# Set instance attributes
		for attrName, value in DEFAULT_SIMULATION_KWARGS.viewitems():
			if attrName in kwargs.viewkeys():
				value = kwargs[attrName]

			setattr(self, "_" + attrName, value)

		unknownKeywords = kwargs.viewkeys() - DEFAULT_SIMULATION_KWARGS.viewkeys()

		if any(unknownKeywords):
			raise SimulationException("Unknown keyword arguments: {}".format(unknownKeywords))

		# Set time variables
		self._simulationStep = 0

		self.randomState = np.random.RandomState(seed = np.uint32(self._seed % np.iinfo(np.uint32).max))

		# divide_cell will fail if _outputDir is no good (e.g. defaulted to
		# None) so catch it *before* running the simulation in case _logToDisk
		# doesn't.
		filepath.makedirs(self._outputDir)

		# Load KB
		sim_data = cPickle.load(open(self._simDataLocation, "rb"))

		# Initialize simulation from fit KB
		self._initialize(sim_data)


	# Link states and processes
	def _initialize(self, sim_data):
		# self._timeStepSec = self._timeStepSec
		self.internal_states = _orderedAbstractionReference(self._internalStateClasses)
		self.external_states = _orderedAbstractionReference(self._externalStateClasses)
		self.processes = _orderedAbstractionReference(self._processClasses)
		self.listeners = _orderedAbstractionReference(self._listenerClasses + DEFAULT_LISTENER_CLASSES)
		self.hooks = _orderedAbstractionReference(self._hookClasses)
		self._initLoggers()
		self._cellCycleComplete = False
		self._isDead = False
		self._finalized = False

		for internal_state in self.internal_states.itervalues():
			internal_state.initialize(self, sim_data)

		for external_state in self.external_states.itervalues():
			external_state.initialize(self, sim_data)

		for process in self.processes.itervalues():
			process.initialize(self, sim_data)

		for listener in self.listeners.itervalues():
			listener.initialize(self, sim_data)

		for hook in self.hooks.itervalues():
			hook.initialize(self, sim_data)

		for internal_state in self.internal_states.itervalues():
			internal_state.allocate()

		for listener in self.listeners.itervalues():
			listener.allocate()

		self._initialConditionsFunction(sim_data)

		self._timeTotal = self.initialTime()

		for hook in self.hooks.itervalues():
			hook.postCalcInitialConditions(self)

		# Make permanent reference to evaluation time listener
		self._evalTime = self.listeners["EvaluationTime"]

		# Perform initial mass calculations
		for state in self.internal_states.itervalues():
			state.calculatePreEvolveStateMass()
			state.calculatePostEvolveStateMass()

		# Update environment state according to the current time in timeseries
		for external_state in self.external_states.itervalues():
			external_state.update()

		# Perform initial listener update
		for listener in self.listeners.itervalues():
			listener.initialUpdate()

		# Start logging
		for logger in self.loggers.itervalues():
			logger.initialize(self)

	def _initLoggers(self):
		self.loggers = collections.OrderedDict()

		if self._logToShell:
			self.loggers["Shell"] = wholecell.loggers.shell.Shell(
				self._shellColumnHeaders
				)

		if self._logToDisk:
			self.loggers["Disk"] = wholecell.loggers.disk.Disk(
				self._outputDir,
				self._overwriteExistingFiles,
				self._logToDiskEvery
				)

	# Run simulation
	def run(self):
		"""
		Run the simulation for the time period specified in `self._lengthSec`
		and then clean up.
		"""

		self.run_incremental(self._lengthSec + self.initialTime())
		self.finalize()

	def run_incremental(self, run_until):
		"""
		Run the simulation for a given amount of time.

		Args:
		    run_until (float): absolute time to run the simulation until. 
		"""

		# Simulate
		while self.time() < run_until and not self._isDead:
			if self._cellCycleComplete:
				self.finalize()
				break

			self._simulationStep += 1

			self._timeTotal += self._timeStepSec

			self._evolveState()

	def finalize(self):
		"""
		Clean up any details once the simulation has finished.
		Specifically, this calls `finalize` in all hooks,
		invokes the simulation's `_divideCellFunction` and then
		shuts down all loggers
		"""

		if not self._finalized:
			# Run post-simulation hooks
			for hook in self.hooks.itervalues():
				hook.finalize(self)

			# Divide mother into daughter cells
			self._divideCellFunction()

			# Finish logging
			for logger in self.loggers.itervalues():
				logger.finalize(self)

			self._finalized = True

	# Calculate temporal evolution
	def _evolveState(self):

		if self._simulationStep <= 1:
			# Update randstreams
			for stateName, state in self.internal_states.iteritems():
				state.seed = self._seedFromName(stateName)
				state.randomState = np.random.RandomState(seed = state.seed)

			for processName, process in self.processes.iteritems():
				process.seed = self._seedFromName(processName)
				process.randomState = np.random.RandomState(seed = process.seed)

		self._adjustTimeStep()

		# Run pre-evolveState hooks
		for hook in self.hooks.itervalues():
			hook.preEvolveState(self)

		# Update queries
		# TODO: context manager/function calls for this logic?
		for i, state in enumerate(self.internal_states.itervalues()):
			t = time.time()
			state.updateQueries()
			self._evalTime.updateQueries_times[i] = time.time() - t

		# Calculate requests
		for i, process in enumerate(self.processes.itervalues()):
			t = time.time()
			process.calculateRequest()
			self._evalTime.calculateRequest_times[i] = time.time() - t

		# Partition states among processes
		for i, state in enumerate(self.internal_states.itervalues()):
			t = time.time()
			state.partition()
			self._evalTime.partition_times[i] = time.time() - t

		# Calculate mass of partitioned molecules
		for state in self.internal_states.itervalues():
			state.calculatePreEvolveStateMass()

		# Update listeners
		for listener in self.listeners.itervalues():
			listener.updatePostRequest()

		# Simulate submodels
		for i, process in enumerate(self.processes.itervalues()):
			t = time.time()
			process.evolveState()
			self._evalTime.evolveState_times[i] = time.time() - t

		# Check that timestep length was short enough
		for process in self.processes.itervalues():
			if not process.wasTimeStepShortEnough():
				raise Exception("The timestep (%.3f) was too long at step %i, failed on process %s" % (self._timeStepSec, self.simulationStep(), str(process.name())))

		# Merge state
		for i, state in enumerate(self.internal_states.itervalues()):
			t = time.time()
			state.merge()
			self._evalTime.merge_times[i] = time.time() - t

		# Calculate mass of partitioned molecules, after evolution
		for state in self.internal_states.itervalues():
			state.calculatePostEvolveStateMass()

		# update environment state
		for state in self.external_states.itervalues():
			state.update()

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
		return np.uint32((self._seed + hash(name)) % np.iinfo(np.uint64).max)
		# return np.uint32((self._seed + self.simulationStep() + hash(name)) % np.iinfo(np.uint64).max)


	def initialTime(self):
		return self._initialTime


	# Save to/load from disk
	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			states = self.internal_states.keys(),
			processes = self.processes.keys()
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStepSec = self.timeStepSec()
			)


	def tableLoad(self, tableReader, tableIndex):
		pass


	def time(self):
		return self._timeTotal


	def simulationStep(self):
		return self._simulationStep


	def timeStepSec(self):
		return self._timeStepSec


	def lengthSec(self):
		return self._lengthSec


	def cellCycleComplete(self):
		self._cellCycleComplete = True


	def _adjustTimeStep(self):
		# Adjust timestep if needed or at a frequency of updateTimeStepFreq regardless
		validTimeSteps = self._maxTimeStep * np.ones(len(self.processes))
		resetTimeStep = False
		for i, process in enumerate(self.processes.itervalues()):
			if not process.isTimeStepShortEnough(self._timeStepSec, self._timeStepSafetyFraction) or self.simulationStep() % self._updateTimeStepFreq == 0:
				validTimeSteps[i] = self._findTimeStep(0., self._maxTimeStep, process.isTimeStepShortEnough)
				resetTimeStep = True
		if resetTimeStep:
			self._timeStepSec = validTimeSteps.min()

	def _findTimeStep(self, minTimeStep, maxTimeStep, checkerFunction):
		N = 10000
		for i in xrange(N):
			candidateTimeStep = minTimeStep + (maxTimeStep - minTimeStep) / 2.
			if checkerFunction(candidateTimeStep, self._timeStepSafetyFraction):
				minTimeStep = candidateTimeStep
				if (maxTimeStep - minTimeStep) / minTimeStep <= 1e-2:
					break
			else:
				maxTimeStep = candidateTimeStep
		if i == N - 1:
			raise Exception, "Timestep adjustment did not converge, last attempt was %f" % (candidateTimeStep)

		return candidateTimeStep

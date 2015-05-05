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
import json
import cPickle
import time

import numpy as np

from wholecell.listeners.evaluation_time import EvaluationTime
from wholecell.containers.bulk_objects_container import BulkObjectsContainer
from wholecell.containers.unique_objects_container import UniqueObjectsContainer

import wholecell.loggers.shell
import wholecell.loggers.disk

from wholecell.io.tablewriter import TableWriter

DEFAULT_SIMULATION_KWARGS = dict(
	seed = 0,
	lengthSec = 3600,
	initialTime = 0,
	logToShell = True,
	logToDisk = False,
	outputDir = None,
	overwriteExistingFiles = False,
	logToDiskEvery = 1,
	kbLocation = None,
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
		"_stateClasses",
		"_processClasses",
		"_initialConditionsFunction",
		)

	# Attributes that may be optionally overwritten by a subclass
	_listenerClasses = ()
	_hookClasses = ()
	_timeStepSec = 1
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
		self.simulationStep = 0

		self.randomState = np.random.RandomState(seed = self._seed)

		# Load KB
		kb = cPickle.load(open(self._kbLocation, "rb"))

		# Initialize simulation from fit KB
		self._initialize(kb)


	# Link states and processes
	def _initialize(self, kb):
		self.states = _orderedAbstractionReference(self._stateClasses)
		self.processes = _orderedAbstractionReference(self._processClasses)
		self.listeners = _orderedAbstractionReference(self._listenerClasses + DEFAULT_LISTENER_CLASSES)
		self.hooks = _orderedAbstractionReference(self._hookClasses)
		self._initLoggers()
		self._cellCycleComplete = False

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

		self._initialConditionsFunction(kb)

		for hook in self.hooks.itervalues():
			hook.postCalcInitialConditions(self)

		# Make permanent reference to evaluation time listener
		self._evalTime = self.listeners["EvaluationTime"]

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
		while self.time() < self._lengthSec + self.initialTime():
			if self._cellCycleComplete:
				break

			self.simulationStep += 1

			self._evolveState()

		# Run post-simulation hooks
		for hook in self.hooks.itervalues():
			hook.finalize(self)

		# Divide mother into daughter cells
		self._divideDaughter()

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

	def _divideDaughter(self):
		# Divide bulk molecules
		molecule_counts = self.states['BulkMolecules'].container.counts()
		d1_bulk_molecules_container = self.states['BulkMolecules'].container.emptyLike()
		d2_bulk_molecules_container = self.states['BulkMolecules'].container.emptyLike()

		d1_counts = self.randomState.binomial(molecule_counts, p = 0.5)
		d2_counts = molecule_counts - d1_counts

		assert all(d1_counts + d2_counts == molecule_counts)

		d1_bulk_molecules_container.countsIs(d1_counts)
		d2_bulk_molecules_container.countsIs(d2_counts)

		# Ensure that special molecules that should not be divded binomially
		# are divided equally between daughter cells
		# TODO: Divide chromosome in some way better!

		# Save data
		os.mkdir(os.path.join(self._outputDir, "Daughter1"))
		os.mkdir(os.path.join(self._outputDir, "Daughter2"))
		d1_table_writer = TableWriter(os.path.join(self._outputDir, "Daughter1", "BulkMolecules"))
		d1_bulk_molecules_container.tableCreate(d1_table_writer)
		d1_bulk_molecules_container.tableAppend(d1_table_writer)

		d2_table_writer = TableWriter(os.path.join(self._outputDir, "Daughter2", "BulkMolecules"))
		d2_bulk_molecules_container.tableCreate(d2_table_writer)
		d2_bulk_molecules_container.tableAppend(d2_table_writer)

		# Divide bulk chromosome
		# chromosome_location_counts = self.states['BulkChromosome'].container.counts()
		# d1_bulk_molecules_container = self.states['BulkChromosome'].container.emptyLike()
		# d2_bulk_molecules_container = self.states['BulkChromosome'].container.emptyLike()
		# import ipdb; ipdb.set_trace()
		# assert all(chromosome_location_counts == 2)

		# d1_counts = chromosome_location_counts / 2
		# d2_counts = chromosome_location_counts / 2

		# assert all(d1_counts + d2_counts == chromosome_location_counts)

		# d1_bulk_chromosome_container.countsIs(d1_counts)
		# d2_bulk_chromosome_container.countsIs(d2_counts)

		# Divide unique molecules
		d1_unique_molecules_container = self.states['UniqueMolecules'].container.emptyLike()
		d2_unique_molecules_container = self.states['UniqueMolecules'].container.emptyLike()

		uniqueMoleculesToDivide = self.states['UniqueMolecules'].uniqueMoleculeDefinitions
		# TODO: We are not dividing dna polymerase!
		uniqueMoleculesToDivide.pop('dnaPolymerase')

		for moleculeName, moleculeAttributeDict in uniqueMoleculesToDivide.iteritems():
			# Get set of molecules to divide and calculate number going to daugher one and daughter two
			moleculeSet = self.states['UniqueMolecules'].container.objectsInCollection(moleculeName)
			n_d1 = self.randomState.binomial(len(moleculeSet), p = 0.5)
			n_d2 = len(moleculeSet) - n_d1
			assert n_d1 + n_d2 == len(moleculeSet)

			# Randomly boolean index molecules in mother such that each daugher gets amount calculated above
			d1_bool = np.zeros(len(moleculeSet), dtype = bool)
			d2_bool = np.zeros(len(moleculeSet), dtype = bool)
			d1_indexes = self.randomState.choice(range(len(moleculeSet)), size = n_d1, replace = False)
			d1_bool[d1_indexes] = True
			d2_bool = np.logical_not(d1_bool)

			d1_dividedAttributesDict = {}
			d2_dividedAttributesDict = {}
			for moleculeAttribute in moleculeAttributeDict.iterkeys():
				d1_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d1_bool]
				d2_dividedAttributesDict[moleculeAttribute] = moleculeSet.attr(moleculeAttribute)[d2_bool]

			d1_unique_molecules_container.objectsNew(moleculeName, n_d1, **d1_dividedAttributesDict)
			d2_unique_molecules_container.objectsNew(moleculeName, n_d2, **d2_dividedAttributesDict)

		# Save data
		d1_table_writer = TableWriter(os.path.join(self._outputDir, "Daughter1", "UniqueMolecules"))
		d1_unique_molecules_container.tableCreate(d1_table_writer)
		d1_unique_molecules_container.tableAppend(d1_table_writer)

		d2_table_writer = TableWriter(os.path.join(self._outputDir, "Daughter2", "UniqueMolecules"))
		d2_unique_molecules_container.tableCreate(d2_table_writer)
		d2_unique_molecules_container.tableAppend(d2_table_writer)

	def _seedFromName(self, name):
		return np.uint32(self._seed + self.simulationStep + hash(name))


	def initialTime(self):
		return self._initialTime


	# Save to/load from disk
	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			states = self.states.keys(),
			processes = self.processes.keys()
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep()
			)


	def tableLoad(self, tableReader, tableIndex):
		pass


	def time(self):
		return self.timeStepSec() * self.simulationStep + self.initialTime()


	def timeStep(self):
		return self.simulationStep


	def timeStepSec(self):
		return self._timeStepSec


	def lengthSec(self):
		return self._lengthSec


	def cellCycleComplete(self):
		self._cellCycleComplete = True

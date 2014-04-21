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

import wholecell.states.mass
import wholecell.states.bulk_molecules
import wholecell.states.unique_molecules
import wholecell.states.chromosome
import wholecell.states.transcripts

STATE_CLASSES = [
	wholecell.states.mass.Mass,
	wholecell.states.bulk_molecules.BulkMolecules,
	wholecell.states.unique_molecules.UniqueMolecules,
	wholecell.states.chromosome.Chromosome,
	wholecell.states.transcripts.Transcripts,
	]

import wholecell.processes.complexation
import wholecell.processes.metabolism
import wholecell.processes.rna_degradation
import wholecell.processes.transcription
import wholecell.processes.translation
import wholecell.processes.free_production
import wholecell.processes.toy_transcription
import wholecell.processes.toy_protein_degradation
import wholecell.processes.toy_replication
import wholecell.processes.replication

PROCESS_CLASSES = [
	wholecell.processes.metabolism.Metabolism,
	wholecell.processes.rna_degradation.RnaDegradation,
	wholecell.processes.transcription.Transcription,
	wholecell.processes.translation.Translation,
	wholecell.processes.free_production.FreeProduction,
	wholecell.processes.toy_transcription.ToyTranscription,
	wholecell.processes.toy_protein_degradation.ToyProteinDegradation,
	wholecell.processes.toy_replication.ToyReplication,
	wholecell.processes.replication.Replication,
	]

STATES = {stateClass.name():stateClass for stateClass in STATE_CLASSES}
PROCESSES = {processClass.name():processClass for processClass in PROCESS_CLASSES}

DEFAULT_STATES = [
	'Mass',
	'BulkMolecules',
	]

DEFAULT_PROCESSES = [
	'Metabolism',
	'RnaDegradation',
	'Transcription',
	'Translation',
	'Replication'
	]

SIM_INIT_ARGS = dict(
	includedStates = None, includedProcesses = None,
	freeMolecules = None,
	lengthSec = None, timeStepSec = None,
	seed = None,
	reconstructKB = False,
	logToShell = True,
	logToDisk = False, outputDir = None, overwriteExistingFiles = False, logToDiskEvery = None
	)

class Simulation(object):
	""" Simulation """

	# Constructors
	def __init__(self, **kwargs):
		import wholecell.utils.rand_stream
		import wholecell.utils.config
		import wholecell.utils.knowledgebase_fixture_manager
		import wholecell.reconstruction.fitter
		import wholecell.loggers.shell
		import wholecell.loggers.disk

		# Make sure the arguments passed are valid
		if (kwargs.viewkeys() - SIM_INIT_ARGS.viewkeys()):
			raise Exception('Unrecognized arguments passed to Simulation.__init__: {}'.format(
					kwargs.viewkeys() - SIM_INIT_ARGS.viewkeys()))

		self._options = SIM_INIT_ARGS.copy()
		self._options.update(kwargs)

		# Set states
		self.includedStates = (self._options['includedStates']
			if self._options['includedStates'] is not None else DEFAULT_STATES)

		# Set processes
		self.includedProcesses = (self._options['includedProcesses']
			if self._options['includedProcesses'] is not None else DEFAULT_PROCESSES)
		
		self.freeMolecules = self._options['freeMolecules']

		if self.freeMolecules is not None:
			self.includedProcesses.append('FreeProduction')

		# Set random seed
		self.seed = self._options['seed']

		# References to simulation objects
		self.randStream = wholecell.utils.rand_stream.RandStream(seed = self.seed)
		self.states = None
		self.processes = None

		# Create KB
		self.kbDir = wholecell.utils.config.SIM_FIXTURE_DIR

		if self._options['reconstructKB'] or not os.path.exists(os.path.join(self.kbDir,'KnowledgeBase.cPickle')):
			kb = wholecell.utils.knowledgebase_fixture_manager.cacheKnowledgeBase(self.kbDir)

		else:
			kb = wholecell.utils.knowledgebase_fixture_manager.loadKnowledgeBase(
				os.path.join(self.kbDir, 'KnowledgeBase.cPickle'))

		# Set time parameters
		self.lengthSec = (self._options['lengthSec']
			if self._options['lengthSec'] is not None else kb.parameters['cellCycleLen'].to('s').magnitude) # Simulation length (s)
		self.timeStepSec = (self._options['timeStepSec']
			if self._options['timeStepSec'] is not None else kb.parameters['timeStep'].to('s').magnitude) # Simulation time step (s)
		self.initialStep = 0
		self.simulationStep = 0

		# Fit KB parameters
		wholecell.reconstruction.fitter.fitKb(kb)
		# TODO: save fit KB and use that instead of saving/loading fit parameters

		# Initialize simulation from fit KB
		self._initialize(kb)

		# Set loggers
		self.loggers = []

		if self._options['logToShell']:
			self.loggers.append(
				wholecell.loggers.shell.Shell()
				)

		if self._options['logToDisk']:
			self.loggers.append(
				wholecell.loggers.disk.Disk(
					self._options['outputDir'],
					self._options['overwriteExistingFiles'],
					self._options['logToDiskEvery']
					)
				)


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
		self._constructStates()
		self._constructProcesses()

		for state in self.states.itervalues():
			state.initialize(self, kb)

		for process in self.processes.itervalues():
			process.initialize(self, kb)

		self._allocateMemory()
		
		from wholecell.reconstruction.initial_conditions import calcInitialConditions

		calcInitialConditions(self, kb)


	# Construct states
	def _constructStates(self):
		self.states = collections.OrderedDict([
			(stateName, STATES[stateName]())
			for stateName in self.includedStates
			])


	# Construct processes
	def _constructProcesses(self):
		self.processes = collections.OrderedDict([
			(processName, PROCESSES[processName]())
			for processName in self.includedProcesses
			])


	# Allocate memory
	def _allocateMemory(self):
		for state in self.states.itervalues():
			state.allocate()


	# -- Run simulation --

	# Run simulation
	def run(self):
		# Calculate initial dependent state
		self._calculateState()

		self._logInitialize()

		while self.time() < self.lengthSec:
			self.simulationStep += 1

			self._evolveState()
			self._logAppend()

		self._logFinalize()


	def _calculateState(self):
		for state in self.states.itervalues():
			state.calculate()


	# Calculate temporal evolution
	def _evolveState(self):
		# Update randstreams
		for stateName, state in self.states.iteritems():
			state.randStream = wholecell.utils.rand_stream.RandStream(
				seed = np.uint64(self.seed + self.simulationStep + hash(stateName))
				)

		for processName, process in self.processes.iteritems():
			process.randStream = wholecell.utils.rand_stream.RandStream(
				seed = np.uint32(self.seed + self.simulationStep + hash(processName))
				)

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

		# Recalculate dependent state
		self._calculateState()


	# --- Logger functions ---
	def _logInitialize(self):
		for logger in self.loggers:
			logger.initialize(self)


	def _logAppend(self):
		for logger in self.loggers:
			logger.append(self)


	def _logFinalize(self):
		for logger in self.loggers:
			logger.finalize(self)

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
	def loadSimulation(cls, stateDir, timePoint, newDir = None, overwriteExistingFiles = False):
		newSim = cls.initFromFile(
			os.path.join(stateDir, 'simOpts.json'),
			logToDisk = newDir is not None,
			overwriteExistingFiles = overwriteExistingFiles,
			outputDir = newDir
			)

		with tables.open_file(os.path.join(stateDir, 'Main.hdf')) as h5file:
			newSim.pytablesLoad(h5file, timePoint)

		for stateName, state in newSim.states.viewitems():
			with tables.open_file(os.path.join(stateDir, stateName + '.hdf')) as h5file:
				state.pytablesLoad(h5file, timePoint)

		newSim.initialStep = timePoint

		# Calculate derived states
		newSim._calculateState() # TODO: add calculate() to State superclass call?

		return newSim


	def options(self):
		return self._options


	def time(self):
		return self.timeStepSec * (self.initialStep + self.simulationStep)


	def timeStep(self):
		return self.initialStep + self.simulationStep

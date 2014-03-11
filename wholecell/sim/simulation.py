#!/usr/bin/env python

"""
Simulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

from __future__ import division

import collections
import tables
import os
import wholecell.utils.config

DEFAULT_PROCESSES = [
	'Complexation',
	'Metabolism',
	'RnaDegradation',
	'Transcription',
	'Translation',
	] # TOKB

SIM_INIT_ARGS = dict(
	includedProcesses = None,
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
		# Make sure the arguments passed are valid
		if (kwargs.viewkeys() - SIM_INIT_ARGS.viewkeys()):
			raise Exception('Unrecognized arguments passed to Simulation.__init__: {}'.format(
					kwargs.viewkeys() - SIM_INIT_ARGS.viewkeys()))

		self._options = SIM_INIT_ARGS.copy()
		self._options.update(kwargs)

		# Set processes
		self.includedProcesses = self._options['includedProcesses'] if self._options['includedProcesses'] is not None else DEFAULT_PROCESSES

		self.freeMolecules = self._options['freeMolecules']

		if self.freeMolecules is not None:
			self.includedProcesses.append('FreeProduction')

		# References to simulation objects
		self.randStream = None
		self.states = None
		self.processes = None

		self._constructRandStream()

		# Set time parameters
		self.lengthSec = self._options['lengthSec'] if self._options['lengthSec'] is not None else 3600. # Simulation length (s) TOKB
		self.timeStepSec = self._options['timeStepSec'] if self._options['timeStepSec'] is not None else 1. # Simulation time step (s) TOKB
		self.initialStep = 0
		self.simulationStep = 0

		# Set random seed
		self.seed = self._options['seed']

		# Create KB
		import wholecell.utils.config
		self.kbDir = wholecell.utils.config.SIM_FIXTURE_DIR

		import cPickle
		if self._options['reconstructKB'] or not os.path.exists(self.kbDir):
			import wholecell.utils.knowledgebase_fixture_manager
			kb = wholecell.utils.knowledgebase_fixture_manager.cacheKnowledgeBase(self.kbDir)
		else:
			import wholecell.utils.knowledgebase_fixture_manager
			kb = wholecell.utils.knowledgebase_fixture_manager.loadKnowledgeBase(
				os.path.join(self.kbDir, 'KnowledgeBase.cPickle'))

		# Fit KB parameters
		import wholecell.reconstruction.fitter
		wholecell.reconstruction.fitter.fitSimulation(kb)
		# TODO: save fit KB and use that instead of saving/loading fit parameters

		# Initialize simulation from fit KB
		self._initialize(kb)

		# Set loggers
		self.loggers = []

		if self._options['logToShell']:
			import wholecell.loggers.shell

			self.loggers.append(
				wholecell.loggers.shell.Shell()
				)

		if self._options['logToDisk']:
			import wholecell.loggers.disk

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
			raise Exception('Caught ValueError; these can be caused by excess commas in the json file, which may not be caught by the syntx checker in your text editor.')

		jsonArgs.update(kwargs)

		return cls(**jsonArgs)


	# Construct random stream
	def _constructRandStream(self):
		import wholecell.utils.rand_stream
		self.randStream = wholecell.utils.rand_stream.RandStream()


	# Link states and processes
	def _initialize(self, kb):
		self._constructStates()
		self._constructProcesses()

		for state in self.states.itervalues():
			state.initialize(self, kb)

		for process in self.processes.itervalues():
			process.initialize(self, kb)

		self._allocateMemory()
		self._calcInitialConditions()


	# Construct states
	def _constructStates(self):
		import wholecell.states.mass
		# import wholecell.states.MetabolicFlux
		import wholecell.states.bulk_molecules
		import wholecell.states.unique_molecules
		import wholecell.states.chromosome
		import wholecell.states.time
		import wholecell.states.rand_stream

		self.states = collections.OrderedDict([
			('Mass',			wholecell.states.mass.Mass()),
			#('MetabolicFlux',	wholecell.sim.state.MetabolicFlux.MetabolicFlux()),
			('BulkMolecules',	wholecell.states.bulk_molecules.BulkMolecules()),
			('UniqueMolecules', wholecell.states.unique_molecules.UniqueMolecules()),
			('Chromosome',		wholecell.states.chromosome.Chromosome()),
			('Time',			wholecell.states.time.Time()),
			('RandStream',		wholecell.states.rand_stream.RandStream())
			])

		self.time = self.states['Time']


	# Construct processes
	def _constructProcesses(self):
		import wholecell.processes.complexation
		import wholecell.processes.metabolism
		import wholecell.processes.rna_degradation
		import wholecell.processes.transcription
		import wholecell.processes.translation
		import wholecell.processes.free_production
		import wholecell.processes.toy_transcription
		import wholecell.processes.toy_protein_degradation

		# TODO: change this so it creates the objects after filtering
		self.processes = collections.OrderedDict([
			('Complexation',		wholecell.processes.complexation.Complexation()),
			('Metabolism',			wholecell.processes.metabolism.Metabolism()),
			('RnaDegradation',		wholecell.processes.rna_degradation.RnaDegradation()),
			('Transcription',		wholecell.processes.transcription.Transcription()),
			('Translation',			wholecell.processes.translation.Translation()),
			('FreeProduction',		wholecell.processes.free_production.FreeProduction()),
			('ToyTranscription',	wholecell.processes.toy_transcription.ToyTranscription()),
			('ToyProteinDegradation',	wholecell.processes.toy_protein_degradation.ToyProteinDegradation())
			])

		# Remove processes not listed as being included
		for process in self.processes.iterkeys():
			if process not in self.includedProcesses:
				self.processes.pop(process)


	# Allocate memory
	def _allocateMemory(self):
		for state in self.states.itervalues():
			state.allocate()


	# Calculate initial conditions
	def _calcInitialConditions(self):
		# Calculate initial conditions
		for state in self.states.itervalues():
			state.calcInitialConditions()


	# -- Run simulation --

	# Run simulation
	def run(self):
		# Calculate initial dependent state
		self._calculateState()

		self._logInitialize()

		while self.time.value < self.lengthSec:
			self.simulationStep += 1

			self._evolveState()
			self._logAppend()

		self._logFinalize()


	def _calculateState(self):
		for state in self.states.itervalues():
			state.calculate()


	# Calculate temporal evolution
	def _evolveState(self):
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
		groupFit = h5file.createGroup(
			h5file.root,
			'fitParameters',
			'Fit parameter values'
			)

		bulkMolecules = self.states['BulkMolecules']
		# TODO: fix all of these/load from fitted KB instead
		# h5file.createArray(groupFit, 'initialDryMass', bulkMolecules.initialDryMass)
		# h5file.createArray(groupFit, 'rnaExp', bulkMolecules.rnaExp)
		# h5file.createArray(groupFit, 'monExp', bulkMolecules.monExp)
		# h5file.createArray(groupFit, 'feistCoreVals', bulkMolecules.feistCoreVals)
		if 'Transcription' in self.processes:
			h5file.createArray(groupFit, 'rnaSynthProb', self.processes['Transcription'].rnaSynthProb)

		groupNames = h5file.createGroup(
			h5file.root,
			'names',
			'State and process names'
			)

		h5file.createArray(groupNames, 'states', [s for s in self.states.viewkeys()])
		
		if self.processes:
			h5file.createArray(groupNames, 'processes', [s for s in self.processes.viewkeys()])

		groupValues = h5file.createGroup(
			h5file.root,
			'values',
			'Non-fit parameter values'
			)

		# h5file.createArray(groupValues, 'molMass', bulkMolecules._molMass)

		# TODO: cache KB


	def pytablesAppend(self, h5file):
		# Included for consistency, eventual features...
		pass


	def pytablesLoad(self, h5file, timePoint):
		group = h5file.get_node('/', 'fitParameters')

		# self.states['BulkMolecules'].initialDryMass = group.initialDryMass.read()
		# self.states['BulkMolecules'].rnaExp[:] = group.rnaExp.read()
		# self.states['BulkMolecules'].monExp[:] = group.monExp.read()
		# self.states['BulkMolecules'].feistCoreVals[:] = group.feistCoreVals.read()
		if 'Transcription' in self.processes:
			self.processes['Transcription'].rnaSynthProb[:] = group.rnaSynthProb.read()

	@classmethod
	def loadSimulation(cls, stateDir, timePoint, newDir = None, overwriteExistingFiles = False):
		newSim = cls.initFromFile(
			os.path.join(stateDir, 'simOpts.json'),
			logToDisk = newDir is not None,
			overwriteExistingFiles = overwriteExistingFiles,
			outputDir = newDir
			)

		with tables.openFile(os.path.join(stateDir, 'Main.hdf')) as h5file:
			newSim.pytablesLoad(h5file, timePoint)

		for stateName, state in newSim.states.viewitems():
			with tables.openFile(os.path.join(stateDir, stateName + '.hdf')) as h5file:
				state.pytablesLoad(h5file, timePoint)

		newSim.initialStep = timePoint

		# Calculate derived states
		newSim._calculateState() # TODO: add calculate() to State superclass call?

		return newSim


	# -- Get, set options and parameters
	def getDynamics(self):
		val = {}
		for state in self.states.itervalues():
			val[state.meta["id"]] = state.getDynamics()
		return val

	@property
	def seed(self):
		return

	@seed.setter
	def seed(self, value):
		if hasattr(self, "randStream"):
			self.randStream.seed = value

	@property
	def randState(self):
		return self.randStream.state

	@randState.setter
	def randState(self, value):
		self.randStream.state = value

	def options(self):
		return self._options

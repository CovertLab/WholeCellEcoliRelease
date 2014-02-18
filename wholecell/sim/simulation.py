#!/usr/bin/env python

"""
Simulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

import numpy
import collections
import tables
import os

DEFAULT_PROCESSES = [
	'Complexation',
	'Metabolism',
	'ProteinMaturation',
	'RnaDegradation',
	'RnaMaturation',
	'Transcription',
	'Translation',
	] # TOKB

KB_PATH = os.path.join('data', 'fixtures', 'KnowledgeBase.cPickle')

# TODO: save/load included processes
# TODO: save/load free molecules

SIM_INIT_ARGS = dict(
	includedProcesses = None,
	freeMolecules = None,
	lengthSec = None, timeStepSec = None,
	seed = None,
	reconstructKB = False, cacheKB = True,
	logToShell = True,
	logToDisk = False, outputDir = None, overwriteExistingFiles = False,
	autoRun = False
	)

class Simulation(object):
	""" Simulation """

	# Constructors
	def __init__(self, **kwargs):
		# Make sure the arguments passed are valid
		if (kwargs.viewkeys() - SIM_INIT_ARGS.viewkeys()):
			raise Exception('Unrecognized arguments passed to Simulation.__init__: {}'.format(
					kwargs.viewkeys() - SIM_INIT_ARGS.viewkeys()))

		options = SIM_INIT_ARGS.copy()
		options.update(kwargs)

		# Set processes
		self.includedProcesses = options['includedProcesses'] if options['includedProcesses'] is not None else DEFAULT_PROCESSES

		self.freeMolecules = options['freeMolecules']

		if self.freeMolecules is not None:
			self.includedProcesses.append('FreeProduction')

		# References to simulation objects
		self.randStream = None
		self.states = None
		self.processes = None

		self.constructRandStream()

		# Set time parameters
		self.lengthSec = options['lengthSec'] if options['lengthSec'] is not None else 3600. # Simulation length (s) TOKB
		self.timeStepSec = options['timeStepSec'] if options['timeStepSec'] is not None else 1. # Simulation time step (s) TOKB
		self.initialStep = 0
		self.simulationStep = 0

		# Set random seed
		self.seed = options['seed']

		# Set KB
		import cPickle
		if options['reconstructKB'] or not os.path.exists(KB_PATH):
			kb = wholecell.reconstruction.knowledgebase.KnowledgeBase(
				dataFileDir = "data/parsed",
				seqFileName = "data/raw/sequence.txt"
				)

			if options['cacheKB']:
				cPickle.dump(kb, open(KB_PATH, "wb"),
					protocol = cPickle.HIGHEST_PROTOCOL)

		else:
			kb = cPickle.load(open(KB_PATH, "rb"))

		self.initialize(kb)

		# Fit model
		import wholecell.utils.fitter
		wholecell.utils.fitter.Fitter.FitSimulation(self, kb)

		# Set loggers
		self.loggers = []

		if options['logToShell']:
			import wholecell.loggers.shell

			self.loggers.append(
				wholecell.loggers.shell.Shell()
				)

		if options['logToDisk']:
			import wholecell.loggers.disk

			self.loggers.append(
				wholecell.loggers.disk.Disk(
					options['outputDir'],
					options['overwriteExistingFiles']
					)
				)

		# Run model (optionally)
		if options['autoRun']:
			self.run()


	@classmethod
	def initFromFile(cls, filePath, **kwargs):
		import json

		try:
			kwargs.update(json.load(open(filePath)))

		except ValueError:
			raise Exception('Caught ValueError; these can be caused by excess commas in the json file, which may not be caught by the syntx checker in your text editor.')


		return cls(**kwargs)


	# Construct random stream
	def constructRandStream(self):
		import wholecell.utils.rand_stream
		self.randStream = wholecell.utils.rand_stream.RandStream()


	# Link states and processes
	def initialize(self, kb):
		self.constructStates()
		self.constructProcesses()

		for state in self.states.itervalues():
			state.initialize(self, kb)

		for process in self.processes.itervalues():
			process.initialize(self, kb)

		self.allocateMemory()
		self.calcInitialConditions()


	# Construct states
	def constructStates(self):
		import wholecell.states.mass
		# import wholecell.states.MetabolicFlux
		import wholecell.states.molecule_counts
		import wholecell.states.time
		import wholecell.states.rand_stream

		self.states = collections.OrderedDict([
			('Mass',			wholecell.states.mass.Mass()),
			#('MetabolicFlux',	wholecell.sim.state.MetabolicFlux.MetabolicFlux()),
			('MoleculeCounts',	wholecell.states.molecule_counts.MoleculeCounts()),
			('Time',			wholecell.states.time.Time()),
			('RandStream',		wholecell.states.rand_stream.RandStream())
			])

		self.time = self.states['Time']


	# Construct processes
	def constructProcesses(self):
		import wholecell.processes.complexation
		import wholecell.processes.metabolism
		import wholecell.processes.protein_maturation
		import wholecell.processes.rna_degradation
		import wholecell.processes.rna_maturation
		import wholecell.processes.transcription
		import wholecell.processes.translation
		import wholecell.processes.free_production

		self.processes = collections.OrderedDict([
			('Complexation',		wholecell.processes.complexation.Complexation()),
			('Metabolism',			wholecell.processes.metabolism.Metabolism()),
			('ProteinMaturation',	wholecell.processes.protein_maturation.ProteinMaturation()),
			('RnaDegradation',		wholecell.processes.rna_degradation.RnaDegradation()),
			('RnaMaturation',		wholecell.processes.rna_maturation.RnaMaturation()),
			('Transcription',		wholecell.processes.transcription.Transcription()),
			('Translation',			wholecell.processes.translation.Translation()),
			('FreeProduction',		wholecell.processes.free_production.FreeProduction())
			])

		# Remove processes not listed as being included
		for process in self.processes.iterkeys():
			if process not in self.includedProcesses:
				self.processes.pop(process)


	# Allocate memory
	def allocateMemory(self):
		for state in self.states.itervalues():
			state.allocate()


	# Calculate initial conditions
	def calcInitialConditions(self):
		# Calculate initial conditions
		for state in self.states.itervalues():
			state.calcInitialConditions()


	# -- Run simulation --

	# Run simulation
	def run(self):
		# Calculate initial dependent state
		self.calculateState()

		self.logInitialize()

		while self.time.value < self.lengthSec:
			self.simulationStep += 1

			self.evolveState()
			self.logAppend()

		self.logFinalize()


	def calculateState(self):
		for state in self.states.itervalues():
			state.calculate()


	# Calculate temporal evolution
	def evolveState(self):
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
		self.calculateState()


	# --- Logger functions ---
	def logInitialize(self):
		for logger in self.loggers:
			logger.initialize(self)


	def logAppend(self):
		for logger in self.loggers:
			logger.append(self)


	def logFinalize(self):
		for logger in self.loggers:
			logger.finalize(self)

	# Save to/load from disk
	def pytablesCreate(self, h5file):
		groupFit = h5file.createGroup(
			h5file.root,
			'fitParameters',
			'Fit parameter values'
			)

		h5file.createArray(groupFit, 'initialDryMass', self.states['MoleculeCounts'].initialDryMass)
		h5file.createArray(groupFit, 'rnaExp', self.states['MoleculeCounts'].rnaExp)
		h5file.createArray(groupFit, 'monExp', self.states['MoleculeCounts'].monExp)
		h5file.createArray(groupFit, 'feistCoreVals', self.states['MoleculeCounts'].feistCoreVals)
		if 'Transcription' in self.processes:
			h5file.createArray(groupFit, 'rnaSynthProb', self.processes['Transcription'].rnaSynthProb)

		groupNames = h5file.createGroup(
			h5file.root,
			'names',
			'State and process names'
			)

		h5file.createArray(groupNames, 'states', [s for s in self.states.viewkeys()])
		h5file.createArray(groupNames, 'processes', [s for s in self.processes.viewkeys()])

		groupValues = h5file.createGroup(
			h5file.root,
			'values',
			'Non-fit parameter values'
			)

		h5file.createArray(groupValues, 'molMass', self.states['MoleculeCounts']._molMass)

		# TODO: save init options


	def pytablesAppend(self, h5file):
		# Included for consistency, eventual features...
		pass


	def pytablesLoad(self, h5file, timePoint):
		group = h5file.get_node('/', 'fitParameters')

		self.states['MoleculeCounts'].initialDryMass = group.initialDryMass.read()
		self.states['MoleculeCounts'].rnaExp[:] = group.rnaExp.read()
		self.states['MoleculeCounts'].monExp[:] = group.monExp.read()
		self.states['MoleculeCounts'].feistCoreVals[:] = group.feistCoreVals.read()
		self.processes['Transcription'].rnaSynthProb[:] = group.rnaSynthProb.read()


	@classmethod
	def loadSimulation(cls, stateDir, timePoint):
		newSim = cls()

		with tables.openFile(os.path.join(stateDir, 'Main.hdf')) as h5file:
			newSim.pytablesLoad(h5file, timePoint)

		for stateName, state in newSim.states.viewitems():
			with tables.openFile(os.path.join(stateDir, stateName + '.hdf')) as h5file:
				try:
					state.pytablesLoad(h5file, timePoint)
				except IndexError:
					raise Exception, 'Time point chosen to load is out of range!\n'

		newSim.initialStep = timePoint

		# Calculate derived states
		newSim.calculateState()

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

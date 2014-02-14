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

# TODO: save/load included processes
# TODO: save/load free molecules

class Simulation(object):
	""" Simulation """

	# Constructor
	def __init__(self, processesToInclude = None, freeMolecules = None):
		self.meta = {
			"options": ["lengthSec", "timeStepSec", "seed"],
			"units": {"lengthSec": "s", "timeStepSec": "s"}
		}

		# Options
		self.lengthSec = 3600						# Simulation length (s) # TOKB
		self.timeStepSec = 1.0						# Simulation time step (s) # TOKB
		if processesToInclude is not None:			# List of processes to include in simulation
			self.processesToInclude = processesToInclude	

		else:
			self.processesToInclude = DEFAULT_PROCESSES

		self.freeMolecules = freeMolecules			# Iterable of tuples describing mol IDs and counts

		if self.freeMolecules is not None:
			self.processesToInclude.append('FreeProduction')

		# Dependent properties
		self.seed = None

		# Rand stream, state, process handles
		self.randStream = None						# Random stream
		self.states = None							# Ordered dict of states
		self.processes = None						# Ordered dict of processes

		self.loggers = []
		self.constructRandStream()

		self.initialStep = 0
		self.simulationStep = 0


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
			if process not in self.processesToInclude:
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
		# Prepare to partition states among processes
		for state in self.states.itervalues():
			state.prepartition()

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
	def loggerAdd(self, logger):
		self.loggers.append(logger)


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
		# TODO: move fitted parameter saving either into the appropriate state,
		# or save the fitted parameters in a knowledge base object

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
	def loadSimulation(cls, kb, stateDir, timePoint):
		newSim = cls()
		newSim.initialize(kb)

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
	def getOptions(self):
		# Initialize output
		val = {"states": {}, "processes": {}}

		# Top-level
		if self.meta.has_key("options"):
			for opt in self.meta["options"]:
				val[opt] = getattr(self, opt)

		# States
		for state in self.states.itervalues():
			val["states"][state.meta["id"]] = state.getOptions()

		# Processes
		for process in self.processes.itervalues():
			val["processes"][process.meta["id"]] = process.getOptions()

		return val

	# Sets options values based on passed in dict
	def setOptions(self, val):
		# Top-level
		opts = list(set(val.keys()).difference(set(["states", "processes"])))
		if len(opts) > 0 and (not self.meta.has_key("options") or not set(opts).issubset(set(self.meta["options"]))):
			invalidOpts = list(set(opts).difference(set(self.meta["options"])))
			raise Exception, "Invalid options:\n -%s" % "".join(invalidOpts)
		for opt in opts:
			setattr(self, opt, val[opt])

		# States
		if val.has_key("states"):
			for key in val["states"].keys():
				state = self.states[key]
				state.setOptions(val["states"][key])

		# Processes
		if val.has_key("processes"):
			for key in val["processes"].keys():
				process = self.processes[key]
				process.setOptions(val["processes"][key])

	def getDynamics(self):
		val = {}
		for state in self.states.itervalues():
			val[state.meta["id"]] = state.getDynamics()
		return val

	@property
	def timeStepSec(self):
		return self._timeStepSec

	@timeStepSec.setter
	def timeStepSec(self, value):
		self._timeStepSec = value
		if hasattr(self, "processes"):
			for process in self.processes.itervalues():
				process.timeStepSec = value

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

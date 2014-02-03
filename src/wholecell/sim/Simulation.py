#!/usr/bin/env python

"""
Simulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

import numpy
import collections

default_processes = ['Complexation', 'Metabolism', 'ProteinMaturation',
					'RnaDegradation', 'RnaMaturation', 'Transcription', 'Translation']

class Simulation(object):
	""" Simulation """

	# Constructor
	def __init__(self, kb, processToInclude = default_processes):
		self.meta = {
			"options": ["lengthSec", "timeStepSec", "seed"],
			"units": {"lengthSec": "s", "timeStepSec": "s"}
		}

		# Options
		self.lengthSec = 50000						# Simulation length (s)
		self.timeStepSec = 1.0						# Simulation time step (s)
		self.processToInclude = processToInclude	# List of processes to include in simulation

		# Dependent properties
		self.seed = None

		# Rand stream, state, process handles
		self.randStream = None						# Random stream
		self.states = None							# Ordered dict of states
		self.time = None							# Time state
		self.processes = None						# Ordered dict of processes

		self.constructRandStream()
		self.constructStates()
		self.constructProcesses()
		self.initialize(kb)
		self.allocateMemory()

	# Construct random stream
	def constructRandStream(self):
		import wholecell.util.randStream
		self.randStream = wholecell.util.randStream.randStream()

	# Construct states
	def constructStates(self):
		import wholecell.sim.state.Mass
		# import wholecell.sim.state.MetabolicFlux
		import wholecell.sim.state.MoleculeCounts
		import wholecell.sim.state.Time
		import wholecell.sim.state.RandStream

		self.states = collections.OrderedDict([
			('Mass',			wholecell.sim.state.Mass.Mass()),
			#('MetabolicFlux',	wholecell.sim.state.MetabolicFlux.MetabolicFlux()),
			('MoleculeCounts',	wholecell.sim.state.MoleculeCounts.MoleculeCounts()),
			('Time',			wholecell.sim.state.Time.Time()),
			('RandStream',		wholecell.sim.state.RandStream.RandStream())
			])

		self.time = self.states["Time"]

	# Construct processes
	def constructProcesses(self):
		import wholecell.sim.process.Complexation
		import wholecell.sim.process.Metabolism
		import wholecell.sim.process.ProteinMaturation
		import wholecell.sim.process.RnaDegradation
		import wholecell.sim.process.RnaMaturation
		import wholecell.sim.process.Transcription
		import wholecell.sim.process.Translation

		self.processes = collections.OrderedDict([
			('Complexation',		wholecell.sim.process.Complexation.Complexation()),
			('Metabolism',			wholecell.sim.process.Metabolism.Metabolism()),
			('ProteinMaturation',	wholecell.sim.process.ProteinMaturation.ProteinMaturation()),
			('RnaDegradation',		wholecell.sim.process.RnaDegradation.RnaDegradation()),
			('RnaMaturation',		wholecell.sim.process.RnaMaturation.RnaMaturation()),
			('Transcription',		wholecell.sim.process.Transcription.Transcription()),
			('Translation',			wholecell.sim.process.Translation.Translation())
			])

	# Link states and processes
	def initialize(self, kb):
		for state in self.states.itervalues():
			state.initialize(self, kb)

		for process in self.processes.itervalues():
			process.initialize(self, kb)

	# Allocate memory
	def allocateMemory(self):
		for state in self.states.itervalues():
			state.allocate()


	# -- Run simulation --

	# Run simulation
	def run(self, loggers = None, timeStep = 0):
		if loggers is None:
			loggers = []

		if timeStep == 0: # stupid hack
			# Calculate initial conditions
			self.calcInitialConditions()

		# Initialize logs
		for logger in loggers:
			logger.initialize(self)

		# Calculate temporal evolution
		for iSec in numpy.arange(self.timeStepSec * (timeStep + 1), self.lengthSec + self.timeStepSec, self.timeStepSec):
			self.time.value = iSec # TODO: update this within Time.calculate()
			self.evolveState()

			# Append logs
			for logger in loggers:
				logger.append(self)

		# Finalize logs
		for logger in loggers:
			logger.finalize(self)

		# print self.randStream.randi(100, 10)

	# Calculate initial conditions
	def calcInitialConditions(self):
		# # Calculate initial conditions
		for state in self.states.itervalues():
			state.calcInitialConditions()

		# Calculate dependent state
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
		for state in self.states.itervalues():
			state.calculate()


	def pytablesCreate(self, h5file):
		# TODO: move fitted parameter saving either into the appropriate state,
		# or save the fitted parameters in a knowledge base object

		group = h5file.createGroup(
			h5file.root,
			'fitParameters',
			'Fit parameter values'
			)

		h5file.createArray(group, 'rnaExp', self.states['MoleculeCounts'].rnaExp)
		h5file.createArray(group, 'monExp', self.states['MoleculeCounts'].monExp)
		h5file.createArray(group, 'feistCoreVals', self.states['MoleculeCounts'].feistCoreVals)
		h5file.createArray(group, 'rnaSynthProb', self.processes['Transcription'].rnaSynthProb)


	def pytablesAppend(self, h5file):
		# Included for consistency, eventual features...
		pass


	def pytablesLoad(self, h5file, timePoint):
		group = h5file.get_node('/', 'fitParameters')

		self.states['MoleculeCounts'].rnaExp[:] = group.rnaExp.read()
		self.states['MoleculeCounts'].monExp[:] = group.monExp.read()
		self.states['MoleculeCounts'].feistCoreVals[:] = group.feistCoreVals.read()
		self.processes['Transcription'].rnaSynthProb[:] = group.rnaSynthProb.read()


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
				state = self.states(key)
				state.setOptions(val["states"][key])

		# Processes
		if val.has_key("processes"):
			for key in val["processes"].keys():
				process = self.processes(key)
				process.setOptions(val["processes"][key])

	# Return parameters as dict
	def getParameters(self):
		# Initialize output
		val = {"states": {}, "processes": {}}

		# States
		for state in self.states.itervalues():
			val["states"][state.meta["id"]] = state.getParameters()

		# Processes
		for process in self.processes.itervalues():
			val["processes"][process.meta["id"]] = process.getParameters()

		return val

	# Sets parameter values based on passed in dict
	def setParameters(self, val):
			# States
			if val.has_key("states"):
				for key in val["states"].keys():
					state = self.states(key)
					state.setParameters(val["states"][key])

			# Processes
			if val.has_key("processes"):
				for key in val["processes"].keys():
					process = self.processes(key)
					process.setOptions(val["processes"][key])

	def getDynamics(self):
		val = {}
		for state in self.states.itervalues():
			val[state.meta["id"]] = state.getDynamics()
		return val

	def setDynamics(self, val):
		for key in val.keys():
			state = self.states(key)
			state.setDynamics(val[key])

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
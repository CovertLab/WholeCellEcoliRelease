#!/usr/bin/env python

"""
Simulation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

import numpy

class Simulation(object):
	""" Simulation """

	# Constructor
	def __init__(self, kb):
		self.meta = {
			"options": ["lengthSec", "timeStepSec", "seed"],
			"units": {"lengthSec": "s", "timeStepSec": "s"}
		}

		# Options
		self.lengthSec = 50000			# Simulation length (s)
		self.timeStepSec = 1.0			# Simulation time step (s)

		# Dependent properties
		self.seed = None

		# Rand stream, state, process handles
		self.randStream = None			# Random stream
		self.states = None				# List of states
		self.time = None				# Time state
		self.processes = None			# List of processes

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
		import wholecell.sim.state.Metabolism
		import wholecell.sim.state.MoleculeCounts
		import wholecell.sim.state.Time

		self.states = [
			wholecell.sim.state.Mass.Mass(),
			wholecell.sim.state.Metabolism.Metabolism(),
			wholecell.sim.state.MoleculeCounts.MoleculeCounts(),
			wholecell.sim.state.Time.Time()
		]

		self.time = self.getState("Time")

	# Construct processes
	def constructProcesses(self):
		import wholecell.sim.process.Complexation
		import wholecell.sim.process.Metabolism
		import wholecell.sim.process.ProteinMaturation
		import wholecell.sim.process.RnaDegradation
		import wholecell.sim.process.RnaMaturation
		import wholecell.sim.process.Transcription
		import wholecell.sim.process.Translation

		self.processes = [
			wholecell.sim.process.Complexation.Complexation(),
			wholecell.sim.process.Metabolism.Metabolism(),
			wholecell.sim.process.ProteinMaturation.ProteinMaturation(),
			wholecell.sim.process.RnaDegradation.RnaDegradation(),
			wholecell.sim.process.RnaMaturation.RnaMaturation(),
			wholecell.sim.process.Transcription.Transcription(),
			wholecell.sim.process.Translation.Translation()
		]

	# Link states and processes
	def initialize(self, kb):
		for state in self.states:
			state.initialize(self, kb)

		for process in self.processes:
			process.initialize(self, kb)

	# Allocate memory
	def allocateMemory(self):
		for state in self.states:
			state.allocate()


	# -- Run simulation --

	# Run simulation
	def run(self, loggers = []):
		# Calculate initial conditions
		self.calcInitialConditions()

		# Initialize logs
		for logger in loggers:
			logger.initialize(self)

		# Calculate temporal evolution
		for iSec in numpy.arange(self.timeStepSec, self.lengthSec + self.timeStepSec, self.timeStepSec):
			self.time.value = iSec
			self.evolveState()

			# Append logs
			for logger in loggers:
				logger.append(self)

		# Finalize logs
		for logger in loggers:
			logger.finalize(self)

	# Calculate initial conditions
	def calcInitialConditions(self):
		# Calculate initial conditions
		self.getState("Time").calcInitialConditions()
		self.getState("MoleculeCounts").calcInitialConditions()
		self.getState("Mass").calculate()
		self.getState("Mass").partition()
		self.getState("Metabolism").calcInitialConditions()

		# Calculate dependent state
		for state in self.states:
			state.calculate()

	# Calculate temporal evolution
	def evolveState(self):
		# Prepare to partition states among processes
		for state in self.states:
			state.prepartition()

		# Partition states among processes
		for state in self.states:
			state.partition()

		# Simulate submodels
		for process in self.processes:
			process.evolveState()

		# Merge state
		for state in self.states:
			state.merge()

		# Recalculate dependent state
		for state in self.states:
			state.calculate()


	# -- Helper methods --

	def getStateIndex(self, stateId):
		for iState in xrange(len(self.states)):
			if self.states[iState].meta["id"] == stateId:
				return iState
		return None

	def getProcessIndex(self, processId):
		for iProcess in xrange(len(self.processes)):
			if self.processes[iProcess].meta["id"] == processId:
				return iProcess
		return None

	def getState(self, stateId):
		iState = self.getStateIndex(stateId)
		if iState != None:
			return self.states[iState]
		return None

	def getProcess(self, processId):
		iProcess = self.getProcessIndex(processId)
		if iProcess != None:
			return self.processes[iProcess]
		return None


	# -- Get, set options and parameters
	def getOptions(self):
		# Initialize output
		val = {"states": {}, "processes": {}}

		# Top-level
		if self.meta.has_key("options"):
			for opt in self.meta["options"]:
				val[opt] = getattr(self, opt)

		# States
		for state in self.states:
			val["states"][state.meta["id"]] = state.getOptions()

		# Processes
		for process in self.processes:
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
				state = self.getState(key)
				state.setOptions(val["states"][key])

		# Processes
		if val.has_key("processes"):
			for key in val["processes"].keys():
				process = self.getProcess(key)
				process.setOptions(val["processes"][key])

	# Return parameters as dict
	def getParameters(self):
		# Initialize output
		val = {"states": {}, "processes": {}}

		# States
		for state in self.states:
			val["states"][state.meta["id"]] = state.getParameters()

		# Processes
		for process in self.processes:
			val["processes"][process.meta["id"]] = process.getParameters()

		return val

	# Sets parameter values based on passed in dict
	def setParameters(self, val):
			# States
			if val.has_key("states"):
				for key in val["states"].keys():
					state = self.getState(key)
					state.setParameters(val["states"][key])

			# Processes
			if val.has_key("processes"):
				for key in val["processes"].keys():
					process = self.getProcess(key)
					process.setOptions(val["processes"][key])

	def getDynamics(self):
		val = {}
		for state in self.states:
			val[state.meta["id"]] = state.getDynamics()
		return val

	def setDynamics(self, val):
		for key in val.keys():
			state = self.getState(key)
			state.setDynamics(val[key])

	@property
	def timeStepSec(self):
		return self._timeStepSec

	@timeStepSec.setter
	def timeStepSec(self, value):
		self._timeStepSec = value
		if hasattr(self, "processes"):
			for process in self.processes:
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
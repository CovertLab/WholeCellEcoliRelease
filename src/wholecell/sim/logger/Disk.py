#!/usr/bin/env python

"""
Disk

Logs whole-cell simulations and metadata to disk.
Also provides a function (load) for reading stored simulation data

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/28/2013
"""

import os
import os.path
import time
import cPickle # TODO: Consider using shelve

import wholecell.sim.logger.Logger

class Disk(wholecell.sim.logger.Logger.Logger):
	""" Disk """

	def __init__(self, metadata = {}, outDir = "", segmentLen = 1000):
		
		if not self.isMetadataValid(metadata):
			raise Exception, "Metadata invalid: %s" % (metadata)

		if outDir == "" or type(outDir) != str:
			raise Exception, "outDir invalid: %s" % (outDir)

		self.metadata = metadata
		self.outDir = outDir
		self.segmentLen = segmentLen

		os.makedirs(self.outDir)

		for f in os.listdir(self.outDir):
			if f.lower().endswith(".cpickle"):
				os.remove(os.path.join(self.outDir, f))

	def initialize(self, sim):
		# -- Metadata --
		self.metadata["startTime"] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
		self.metadata["endTime"] = []
		self.metadata["lengthSec"] = []
		self.metadata["timeStepSec"] = []
		self.metadata["segmentLen"] = []

		# -- Save initial state --
		# Setup segment, step counters
		self.iSegment = 0
		self.iStep = 1

		# Allocate memory
		self.allocateMemory(sim, 1)

		# Copy data from states
		self.copyDataFromStates(sim)

		# Save initial state to disk
		self.saveSegmentToDisk()

		# -- Setup to save dynamics --
		self.iSegment = 1
		self.iStep = 0

		# Allocate memory
		self.allocateMemory(sim, self.segmentLen)

	def append(self, sim):
		# Increment step counter
		self.iStep += 1

		# Copy data from states
		self.copyDataFromStates(sim)

		# Rotate segment
		if self.iStep == self.segmentLen:
			# Save segment to disk
			self.saveSegmentToDisk()

			# Increment segment counter
			self.iStep = 0
			self.iSegment += 1

	def finalize(self, sim):
		# -- Save final segment --
		if self.iStep > 0:
			# Contract segment
			self.contractSegment(sim)

			# Save last segment to disk
			self.saveSegmentToDisk

		# -- Metadata --
		# Record
		self.metadata["endTime"] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
		self.metadata["lengthSec"] = sim.getState("Time").value
		self.metadata["timeStepSec"] = sim.timeStepSec

		# Save
		self.saveMetadata(sim.getOptions(), sim.getParameters())

		# -- Clear log --
		del self.iSegment
		del self.iStep
		del self.stateLog
		del self.randStreamLog


	# TODO: Implement. Check for the following:
	# - Name, description
	# - Investigator
	# - Revision
	# - Username, hostname, ip
	def isMetadataValid(self, metadata):
		return True

	def allocateMemory(self, sim, nSteps):
		# State
		self.stateLog = {}

		for state in sim.states:
			stateId = state.meta["id"]
			self.stateLog[stateId] = {}
			for prop in s.meta["dynamics"]:
				self.stateLog[stateId][prop] = numpy.zeros(numpy.shape(state[prop]) + (nSteps,))

		# Random stream
		self.randStreamLog = []

	def contractSegment(self, sim):
		# State
		for state in sim.states:
			for prop in state.meta["dynamics"]:
				self.stateLog[stateId][prop] = self.stateLog[stateId][prop][:, :, 0:(self.iStep + 1)]

	def copyDataFromStates(self, sim):
		# State
		for state in sim.states:
			stateId = state.meta["id"]
			for prop in state.meta["dynamics"]:
				self.stateLog[stateId][prop][:, :, self.iStep] = getattr(state, prop)

		# Random stream
		self.randStreamLog.append(sim.randStream.state)

	def saveSegmentToDisk(self):
		# Dynamics
		cPickle.dump(self.stateLog, os.path.join(self.outDir, "state-%d.cPickle" % (self.iSegment)), protocol = cPickle.HIGHEST_PROTOCOL)

		# Random stream
		cPickle.dump(self.randStreamLog, os.path.join(self.outDir, "randStream-%d.cPickle" % (self.iSegment)), protocol = cPickle.HIGHEST_PROTOCOL)

	def saveMetadata(self, options, parameters):
		metadata = self.metadata
		cPickle.dump(metadata, os.path.join(self.outDir, "metadata.cPickle"), protocol = cPickle.HIGHEST_PROTOCOL)
		cPickle.dump(options, os.path.join(self.outDir, "options.cPickle"), protocol = cPickle.HIGHEST_PROTOCOL)
		cPickle.dump(parameters, os.path.join(self.outDir, "parameters.cPickle"), protocol = cPickle.HIGHEST_PROTOCOL)

	# TODO: Write this method after we actually have some data to work with
	@classmethod
	def load(cls, dir, state, prop):
		pass
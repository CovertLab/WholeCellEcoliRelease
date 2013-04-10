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
import numpy

import wholecell.sim.logger.Logger

class Disk(wholecell.sim.logger.Logger.Logger):
	""" Disk """

	def __init__(self, metadata = None, outDir = None, segmentLen = 1000):
		
		if not self.isMetadataValid(metadata):
			raise Exception, "Metadata invalid: %s" % (metadata)

		if outDir == None or type(outDir) != str:
			raise Exception, "outDir invalid: %s" % (outDir)

		self.metadata = metadata or {}
		self.outDir = outDir
		self.segmentLen = segmentLen

		if not os.path.exists(self.outDir):
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
		self.metadata["segmentLen"] = self.segmentLen

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
			for prop in state.meta["dynamics"]:
				self.stateLog[stateId][prop] = [] #numpy.zeros(numpy.shape(getattr(state, prop)) + (nSteps,))

		# Random stream
		self.randStreamLog = []

	# TODO: Is this function necessary since we're just doing dstack?
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
				if self.iStep == 1:
					self.stateLog[stateId][prop] = getattr(state, prop)
				else:
					self.stateLog[stateId][prop] = numpy.dstack((self.stateLog[stateId][prop], getattr(state, prop)))

		# Random stream
		self.randStreamLog.append(sim.randStream.state)

	def saveSegmentToDisk(self):
		# Dynamics
		cPickle.dump(self.stateLog, open(os.path.join(self.outDir, "state-%d.cPickle" % (self.iSegment)), "w"), protocol = cPickle.HIGHEST_PROTOCOL)

		# Random stream
		cPickle.dump(self.randStreamLog, open(os.path.join(self.outDir, "randStream-%d.cPickle" % (self.iSegment)), "w"), protocol = cPickle.HIGHEST_PROTOCOL)

	def saveMetadata(self, options, parameters):
		metadata = self.metadata
		cPickle.dump(metadata, open(os.path.join(self.outDir, "metadata.cPickle"), "w"), protocol = cPickle.HIGHEST_PROTOCOL)
		cPickle.dump(options, open(os.path.join(self.outDir, "options.cPickle"), "w"), protocol = cPickle.HIGHEST_PROTOCOL)
		cPickle.dump(parameters, open(os.path.join(self.outDir, "parameters.cPickle"), "w"), protocol = cPickle.HIGHEST_PROTOCOL)

	@classmethod
	def load(cls, folder, state, prop):
		# Load metadata
		md = cPickle.load(open(os.path.join(folder, "metadata.cPickle"), "r"))

		# Load first time point
		tmp = cPickle.load(open(os.path.join(folder, "state-0.cPickle"), "r"))

		# Store initial segment
		value = tmp[state][prop]

		# Load and store subsequent segments
		for iTime in xrange(int(md["timeStepSec"]), int(md["lengthSec"]), int(md["segmentLen"] * md["timeStepSec"])):
			iSegment = (iTime - md["timeStepSec"]) / (md["timeStepSec"] * md["segmentLen"]) + 1
			tmp = cPickle.load(open(os.path.join(folder, "state-%d.cPickle" % iSegment), "r"))
			value = numpy.dstack((value, tmp[state][prop]))

		return value
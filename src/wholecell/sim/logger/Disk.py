#!/usr/bin/env python

"""
Disk

Logs whole-cell simulations and metadata to disk using pytables.
Also provides a function (load) for reading stored simulation data

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/3/2013
"""

import os
import time
import tables

import wholecell.sim.logger.Logger

# TODO: let loaded simulation resume logging in a copied file

class Disk(wholecell.sim.logger.Logger.Logger):
	""" Disk """

	def __init__(self, outDir = None):
		self.outDir = outDir
		self.stateFiles = {}
		self.mainFile = None


	def initialize(self, sim):
		if self.outDir is None:
			self.outDir = os.path.join('out', self.currentTimeAsDir())
		
		os.makedirs(self.outDir)

		self.mainFile = tables.open_file(
			os.path.join(self.outDir, 'Main.hdf'),
			mode = "w",
			title = "Main simulation file"
			)
		
		# Metadata
		self.mainFile.root._v_attrs.startTime = self.currentTimeAsString()
		self.mainFile.root._v_attrs.timeStepSec = sim.timeStepSec

		# Create tables
		self.createTables(sim)

		# Save initial state
		self.copyDataFromStates(sim)


	def append(self, sim):
		self.copyDataFromStates(sim)


	def finalize(self, sim):
		# Metadata
		self.mainFile.root._v_attrs.lengthSec = sim.states['Time'].value
		self.mainFile.root._v_attrs.endTime = self.currentTimeAsString()

		# Close file
		self.mainFile.close()

		for stateFile in self.stateFiles.viewvalues():
			stateFile.close()


	def createTables(self, sim):
		sim.pytablesCreate(self.mainFile)

		for stateName, state in sim.states.viewitems():
			stateFile = tables.open_file(
				os.path.join(self.outDir, stateName + '.hdf'),
				mode = "w",
				title = stateName + " state file"
				)

			state.pytablesCreate(stateFile)

			self.stateFiles[state] = stateFile


	def copyDataFromStates(self, sim):
		sim.pytablesAppend(self.mainFile)

		for state, stateFile in self.stateFiles.viewitems():
			state.pytablesAppend(stateFile)


	@staticmethod
	def currentTimeAsString():
		return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


	@staticmethod
	def currentTimeAsDir():
		# Variant timestamp format that should be valid across OSes
		return time.strftime("sim%Y.%m.%d-%H.%M.%S", time.localtime())

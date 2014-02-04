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

	def __init__(self, outDir):
		self.outDir = outDir
		self.h5file = None

		if not os.path.exists(outDir):
			os.makedirs(outDir)


	def initialize(self, sim):
		self.h5file = tables.open_file(
			os.path.join(self.outDir, 'state.hdf'),
			mode = "w",
			title = "Single simulation"
			)
		
		# Metadata
		self.h5file.root._v_attrs.startTime = self.currentTimeAsString()
		self.h5file.root._v_attrs.timeStepSec = sim.timeStepSec
		# self.h5file.root._v_attrs.options = sim.getOptions()
		# self.h5file.root._v_attrs.parameters = sim.getParameters()

		# Create tables
		self.createTables(sim)

		# Save initial state
		self.copyDataFromStates(sim)


	def append(self, sim):
		self.copyDataFromStates(sim)


	def finalize(self, sim):
		# Metadata
		self.h5file.root._v_attrs.lengthSec = sim.states['Time'].value
		self.h5file.root._v_attrs.endTime = self.currentTimeAsString()

		# Close file
		self.h5file.close()


	def createTables(self, sim):
		sim.pytablesCreate(self.h5file)

		for state in sim.states.itervalues():
			state.pytablesCreate(self.h5file)


	def copyDataFromStates(self, sim):
		sim.pytablesAppend(self.h5file)

		for state in sim.states.itervalues():
			state.pytablesAppend(self.h5file)


	@staticmethod
	def currentTimeAsString():
		time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

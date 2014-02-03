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
import cPickle

import wholecell.sim.logger.Logger
import wholecell.util.PickleHelper

# TODO: add support for all States
# TODO: add loading code
# TODO: pickle the initial simulation

class Disk(wholecell.sim.logger.Logger.Logger):
	""" Disk """

	def __init__(self, outDir = None):
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

		# Pickle initial simulation
		wholecell.util.PickleHelper.registerInstanceMethods()
		cPickle.dump(
			sim,
			open(os.path.join(self.outDir, 'simulation.cPickle'), 'wb'),
			protocol = cPickle.HIGHEST_PROTOCOL
			)


	def append(self, sim):
		self.copyDataFromStates(sim)


	def finalize(self, sim):
		# Metadata
		self.h5file.root._v_attrs.lengthSec = sim.getState('Time').value
		self.h5file.root._v_attrs.endTime = self.currentTimeAsString()

		# Close file
		self.h5file.close()


	def createTables(self, sim):
		for state in sim.states:
			state.pytablesCreate(self.h5file)


	def copyDataFromStates(self, sim):
		for state in sim.states:
			state.pytablesAppend(self.h5file)


	@staticmethod
	def currentTimeAsString():
		time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

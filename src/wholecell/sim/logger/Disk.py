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
import os.path
import time
import tables

import wholecell.sim.logger.Logger

# TODO: add support for all States
# TODO: add loading code
# TODO: pickle the initial simulation
# TODO: remove the "sim" reference in copy/create?

class Disk(wholecell.sim.logger.Logger.Logger):
	""" Disk """

	def __init__(self, outFile = None, overwrite = False):
		self.outFile = outFile
		self.h5file = None

		if os.path.exists(self.outFile) and not overwrite:
			raise Exception, "File exists. To overwrite, specify overwrite = True."

		else:
			self.h5file = tables.open_file(outFile, mode = "w", title = "Single simulation")


	def initialize(self, sim):
		# -- Metadata --
		self.h5file.root._v_attrs.startTime = self.currentTimeAsString()
		self.h5file.root._v_attrs.timeStepSec = sim.timeStepSec
		# self.h5file.root._v_attrs.options = sim.getOptions()
		# self.h5file.root._v_attrs.parameters = sim.getParameters()

		# -- Create tables --
		self.createTables(sim)

		# -- Save initial state --
		self.copyDataFromStates(sim)


	def append(self, sim):
		self.copyDataFromStates(sim)


	def finalize(self, sim):
		# -- Metadata --
		self.h5file.root._v_attrs.lengthSec = sim.getState('Time').value
		self.h5file.root._v_attrs.endTime = self.currentTimeAsString()

		# -- Close file --
		self.h5file.close()


	def createTables(self, sim):
		for state in sim.states:
			state.pytablesCreate(self.h5file, sim)


	def copyDataFromStates(self, sim):
		# State
		for state in sim.states:
			state.pytablesAppend(self.h5file, sim)


	@staticmethod
	def currentTimeAsString():
		time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

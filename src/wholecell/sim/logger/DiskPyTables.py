#!/usr/bin/env python

"""
DiskPyTables

Logs whole-cell simulations and metadata to disk using pytables.
Also provides a function (load) for reading stored simulation data

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/3/2013
"""

import os
import os.path
import time
import tables

import wholecell.sim.logger.Logger

class DiskPyTables(wholecell.sim.logger.Logger.Logger):
	""" DiskPyTables """

	def __init__(self, metadata = None, outFile = None, overwrite = False):
		
		if not self.isMetadataValid(metadata):
			raise Exception, "Metadata invalid: %s" % (metadata)

		if outFile == None or type(outFile) != str:
			raise Exception, "outFile invalid: %s" % (outFile)

		self.metadata = metadata or {}
		self.outFile = outFile
		self.h5file = None

		if os.path.exists(self.outFile) and not overwrite:
			raise Exception, "File exists. To overwrite, specify overwrite = True."
		else:
			self.h5file = tables.open_file(outFile, mode = "w", title = "Single simulation")

	def initialize(self, sim):
		# -- Metadata --

		self.metadata["startTime"] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
		self.metadata["endTime"] = []
		self.metadata["lengthSec"] = []
		self.metadata["timeStepSec"] = []

		# -- Create tables --
		self.createTables(sim)

		# -- Save initial state --		
		self.copyDataFromStates(sim)

	def append(self, sim):
		self.copyDataFromStates(sim)

	def finalize(self, sim):

		# -- Metadata --
		# Record
		self.metadata["endTime"] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
		self.metadata["lengthSec"] = sim.getState("Time").value
		self.metadata["timeStepSec"] = sim.timeStepSec

		# Save
		self.saveMetadata(sim.getOptions(), sim.getParameters())

		# -- Close file --
		self.h5file.close()


	# TODO: Implement. Check for the following:
	# - Name, description
	# - Investigator
	# - Revision
	# - Username, hostname, ip
	def isMetadataValid(self, metadata):
		return True

	def createTables(self, sim):
		for state in sim.states:
			if hasattr(state, "pytablesCreate"):
				state.pytablesCreate(self.h5file, sim)

	def copyDataFromStates(self, sim):
		# State
		for state in sim.states:
			if hasattr(state, "pytablesAppend"):
				state.pytablesAppend(self.h5file, sim)

	def saveMetadata(self, options, parameters):
		metadata = self.metadata

		self.h5file.root._v_attrs.metadata = metadata
		self.h5file.root._v_attrs.options = options
		self.h5file.root._v_attrs.parameters = parameters

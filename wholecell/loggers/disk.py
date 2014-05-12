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

from __future__ import division

import os
import time
import json
import shutil
import itertools

import tables

import wholecell.loggers.logger

# TODO: let loaded simulation resume logging in a copied file

DEFAULT_LOG_FREQUENCY = 1

class Disk(wholecell.loggers.logger.Logger):
	""" Disk """

	def __init__(self, outDir = None, allowOverwrite = False, logEvery = None):
		self.outDir = outDir
		self.allowOverwrite = allowOverwrite
		self.logEvery = logEvery if logEvery is not None else DEFAULT_LOG_FREQUENCY

		self.saveFiles = {}
		self.mainFile = None
		self.logStep = None


	def initialize(self, sim):
		self.logStep = 0

		if self.outDir is None:
			self.outDir = os.path.join('out', self.currentTimeAsDir())
		
		try:
			os.makedirs(self.outDir)

		except OSError as e:
			if self.allowOverwrite:
				pass

			else:
				raise

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
		self.copyData(sim)

		# Save simulation definition
		json.dump(
			sim.options().toDict(),
			open(os.path.join(self.outDir, 'simOpts.json'), 'w'),
			sort_keys = True, indent=4, separators=(',', ': ')
			)

		# Save KB
		# shutil.copy(sim.kbPath, self.outDir)
		# TODO: reinstate


	def createTables(self, sim):
		expectedRows = int(sim.lengthSec/sim.timeStepSec)

		sim.pytablesCreate(self.mainFile, expectedRows)

		for name, obj in itertools.chain(sim.states.viewitems(), sim.listeners.viewitems()):
			saveFile = tables.open_file(
				os.path.join(self.outDir, name + '.hdf'),
				mode = "w",
				title = name + " simulation data file"
				)

			obj.pytablesCreate(saveFile, expectedRows)

			self.saveFiles[obj] = saveFile


	def append(self, sim):
		self.logStep += 1

		if self.logStep % self.logEvery == 0:
			self.copyData(sim)


	def finalize(self, sim):
		# Metadata
		self.mainFile.root._v_attrs.lengthSec = sim.time()
		self.mainFile.root._v_attrs.endTime = self.currentTimeAsString()

		# Close file
		self.mainFile.close()

		for saveFile in self.saveFiles.viewvalues():
			saveFile.close()


	def copyData(self, sim):
		sim.pytablesAppend(self.mainFile)

		for obj, saveFile in self.saveFiles.viewitems():
			obj.pytablesAppend(saveFile)


	@staticmethod
	def currentTimeAsString():
		return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


	@staticmethod
	def currentTimeAsDir():
		# Variant timestamp format that should be valid across OSes
		return time.strftime("sim%Y.%m.%d-%H.%M.%S", time.localtime())

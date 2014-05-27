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
from numpy import log10

import wholecell.loggers.logger
from wholecell.utils.constants import OUTPUT_DIRECTORY

# TODO: let loaded simulation resume logging in a copied file

DEFAULT_LOG_FREQUENCY = 1

class Disk(wholecell.loggers.logger.Logger):
	""" Disk """

	def __init__(self, outDir = None, allowOverwrite = False, logEvery = None):
		if outDir is None:
			self.outDir = timestampedOutputDirectoryNames()[0]

		else:
			self.outDir = outDir
		
		try:
			os.makedirs(self.outDir)

		except OSError as e:
			if self.allowOverwrite:
				pass

			else:
				raise

		self.allowOverwrite = allowOverwrite
		self.logEvery = logEvery if logEvery is not None else DEFAULT_LOG_FREQUENCY

		self.saveFiles = {}
		self.mainFile = None
		self.logStep = 0


	def initialize(self, sim):
		self.mainFile = tables.open_file(
			os.path.join(self.outDir, 'Main.hdf'),
			mode = "w",
			title = "Main simulation file"
			)
		
		# Metadata
		self.mainFile.root._v_attrs.startTime = currentTimeAsString()
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
		self.mainFile.root._v_attrs.endTime = currentTimeAsString()

		# Close file
		self.mainFile.close()

		for saveFile in self.saveFiles.viewvalues():
			saveFile.close()


	def copyData(self, sim):
		sim.pytablesAppend(self.mainFile)

		for obj, saveFile in self.saveFiles.viewitems():
			obj.pytablesAppend(saveFile)


def currentTimeAsString():
	return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())



def currentTimeAsDir():
	# Variant timestamp format that should be valid across OSes
	return time.strftime("sim%Y.%m.%d-%H.%M.%S", time.localtime())


def timestampedOutputDirectoryNames(n = 1):
	mainPath = os.path.join(OUTPUT_DIRECTORY, currentTimeAsDir())

	nDigits = int(log10(n)) + 1

	return [
		os.path.join(mainPath, str(simulationIndex).zfill(nDigits))
		for simulationIndex in xrange(n)
		]

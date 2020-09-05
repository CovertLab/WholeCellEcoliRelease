"""
Disk

Logs whole-cell simulations and metadata to disk.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/3/2013
"""

from __future__ import absolute_import, division, print_function

import os
import time
import itertools

import wholecell.loggers.logger
from wholecell.io.tablewriter import TableWriter
import six

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
		self.logStep = 0


	def initialize(self, sim):
		self.mainFile = TableWriter(os.path.join(self.outDir, "Main"))

		# Metadata
		self.mainFile.writeAttributes(
			initialTime = sim.initialTime(),
			startTime = currentTimeAsString(),
			)

		# Create tables
		self.createTables(sim)

		# Save initial state
		self.copyData(sim)

		# TODO: save simulation settings


	def createTables(self, sim):
		sim.tableCreate(self.mainFile)

		# TODO: separate checkpointing and logging
		for name, obj in itertools.chain(six.viewitems(sim.internal_states), six.viewitems(sim.external_states), six.viewitems(sim.listeners)):
			saveFile = TableWriter(os.path.join(self.outDir, name))

			obj.tableCreate(saveFile)

			self.saveFiles[obj] = saveFile


	def append(self, sim):
		self.logStep += 1

		if self.logStep % self.logEvery == 0:
			self.copyData(sim)


	def finalize(self, sim):
		# Metadata
		self.mainFile.writeAttributes(
			lengthSec = sim.time(),
			endTime = currentTimeAsString()
			)

		# Close files
		self.mainFile.close()

		for saveFile in six.viewvalues(self.saveFiles):
			saveFile.close()


	def copyData(self, sim):
		sim.tableAppend(self.mainFile)

		for obj, saveFile in six.viewitems(self.saveFiles):
			obj.tableAppend(saveFile)


def currentTimeAsString():
	return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

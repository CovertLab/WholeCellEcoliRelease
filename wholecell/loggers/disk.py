"""
Disk

Logs whole-cell simulations and metadata to disk.
"""

from __future__ import annotations

import itertools
import os
import time
from typing import Any, Dict

import wholecell.loggers.logger
from wholecell.io.tablewriter import TableWriter
import six

# TODO: let loaded simulation resume logging in a copied file

DEFAULT_LOG_FREQUENCY = 1

class Disk(wholecell.loggers.logger.Logger):
	"""
	Writes state and listener data to disk.  Responsible for creating and
	appending to data tables as well as writing attributes.
	"""

	def __init__(self, outDir: str, logEvery: int = DEFAULT_LOG_FREQUENCY):
		"""
		Args:
			outDir: directory path to save data to
			logEvery: how frequently (# of timesteps) data is saved to disk
				(ie 2 will write data for every other time step)
		"""

		self.outDir = outDir
		self.logEvery = logEvery

		self.saveFiles = {}  # type: Dict[Any, TableWriter]
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

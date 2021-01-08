#!/usr/bin/env python

"""
Shell

Prints a very brief summary of a whole-cell simulation to standard output
"""

from __future__ import annotations

import datetime
import numpy as np
import os
import sys
from typing import Iterable, Optional

from six.moves import range, zip

import wholecell.loggers.logger
from wholecell.utils.py3 import monotonic_seconds


SPACER = "  "

class Shell(wholecell.loggers.logger.Logger):
	"""
	Displays a simple summary of the simulation state to the shell as the
	simulation progresses.  Optionally saves the output in a log file.
	Logged values are added to columnSpecs by calling registerLoggedQuantity()
	in a Listener class.
	"""

	def __init__(self,
			columnHeaders: Iterable[str],
			output_dir: Optional[str] = None,
			):
		"""
		Args:
			columnHeaders: header IDs for values added to columnSpecs to
				display at each time step
			output_dir: if not None, will also save the output to a log file
				in this directory in addition to logging to the shell
		"""

		self.iterFreq = 1
		self.headerFreq = 50

		# Can also be populated by Listener classes calling registerLoggedQuantity()
		self.columnSpecs = [
			{"header": "Time (s)", "target": "Simulation", "property": "time", "length": 8, "format": ".2f", "sum": False},
			]

		self.columnHeaders = columnHeaders

		self.columns = None
		self._header = None
		self._headerBoundary = None

		self.nLines = -1
		self.startTime = monotonic_seconds()
		self.startSimTime = 0

		self.log_file = None
		if output_dir:
			self.log_file = open(os.path.join(output_dir, 'shell.log'), 'w')


	def initialize(self, sim):
		self.columns = []

		for header in self.columnHeaders:
			for spec in self.columnSpecs:
				if spec["header"].replace("\n", " ") == header:
					self.columns.append(spec)
					break

			else:
				raise Exception("Could not find column named {}".format(header))

		# Build the header
		self._buildHeader()

		# Collect Metadata
		self.nLines = -1
		self.startTime = monotonic_seconds()
		self.startSimTime = sim.time()

		# Print initial state
		self.append(sim)


	def printHeaders(self):
		if self.nLines > 0:
			self.write(self._headerBoundary)

		self.write(self._header)

		self.write(self._headerBoundary)


	def _buildHeader(self):
		columnHeaders = [columnSpec["header"] for columnSpec in self.columns]
		cellSizes = [columnSpec["length"] for columnSpec in self.columns]

		# Break the headers at newline characters
		columnHeaderLines = [
			columnHeader.splitlines() for columnHeader in columnHeaders
			]

		# Update the cell size to be at least the header width
		cellSizes = [
			max(cellSize, max(len(line) for line in lines))
			for cellSize, lines in zip(cellSizes, columnHeaderLines)
			]

		# Rearrange the header lines
		headerLines = []
		for headerLineIndex in range(max(len(lines) for lines in columnHeaderLines)):
			line = []
			for lines in columnHeaderLines:
				if len(lines) > headerLineIndex:
					line.append(lines[headerLineIndex])

				else:
					line.append("")

			headerLines.append(line)

		# Build the header

		strings = []

		for headers in headerLines:
			string = []
			for columnIndex, (columnSize, columnHeader) in enumerate(zip(cellSizes, headers)):
				if columnIndex > 0:
					string.append(SPACER)

				string.append(("%" + str(columnSize) + "s") % columnHeader)

			strings.append(''.join(string))

		self._header = '\n'.join(strings) + '\n'

		string = []
		for columnIndex, columnSize in enumerate(cellSizes):
			if columnIndex > 0:
				string.append(SPACER)

			string.append(("%" + str(columnSize) + "s") % ("=" * columnSize))

		self._headerBoundary = ''.join(string) + '\n'

		# Update cell sizes

		for columnSpec, cellSize in zip(self.columns, cellSizes):
			columnSpec["length"] = cellSize


	def append(self, sim):
		self.nLines += 1

		if self.nLines % self.iterFreq != 0:
			return

		if self.nLines % self.headerFreq == 0:
			self.printHeaders()

		for iColumn in range(len(self.columns)):
			column = self.columns[iColumn]

			if iColumn > 0:
				self.write(SPACER)

			if column["target"] == "Simulation":
				target = sim

			else:
				targetType, targetName = column["target"].split(":")

				target = {
					"State":sim.internal_states,
					"Process":sim.processes,
					"Listener":sim.listeners
					}[targetType][targetName]

			value = getattr(target, column["property"])

			if callable(value):
				value = value()

			if column["sum"]:
				value = np.sum(value)

			self.write(("%" + str(column["length"]) + column["format"]) % value)

		self.write("\n")


	def finalize(self, sim):
		# Print summary
		self.write("\n")
		self.write("Simulation finished:\n")

		simTime = sim.time()
		simLength = simTime - self.startSimTime
		runtime = monotonic_seconds() - self.startTime

		self.write(" - Sim length: {}\n".format(hms(simLength)))
		self.write(" - Sim end time: {}\n".format(hms(simTime)))
		self.write(" - Runtime: {}\n".format(hms(runtime)))

		self.write("\n")

	def write(self, text):
		sys.stdout.write(text)
		if self.log_file:
			self.log_file.write(text)

def hms(seconds):
	"""Format a time interval of seconds as [days] h:mm:ss."""
	# Rounding gets e.g. '0:08:29' instead of '0:08:28.809659'.
	delta = datetime.timedelta(seconds=round(seconds))
	return str(delta)

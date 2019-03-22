#!/usr/bin/env python

"""
Shell

Prints a very brief summary of a whole-cell simulation to standard output

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

from __future__ import division

import time
import sys
import numpy as np
from itertools import izip

import wholecell.loggers.logger

SPACER = "  "

class Shell(wholecell.loggers.logger.Logger):
	""" Shell """

	def __init__(self, columnHeaders):
		self.iterFreq = 1
		self.headerFreq = 50

		self.columnSpecs = [
			{"header": "Time (s)", "target": "Simulation", "property": "time", "length": 8, "format": ".2f", "sum": False},
			]

		self.columnHeaders = columnHeaders

		self.columns = None
		self._header = None
		self._headerBoundary = None


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
		self.startTime = time.time()

		# Print initial state
		self.append(sim)


	def printHeaders(self):
		if self.nLines > 0:
			sys.stdout.write(self._headerBoundary)

		sys.stdout.write(self._header)

		sys.stdout.write(self._headerBoundary)


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
			for cellSize, lines in izip(cellSizes, columnHeaderLines)
			]

		# Rearrange the header lines
		headerLines = []
		for headerLineIndex in xrange(max(len(lines) for lines in columnHeaderLines)):
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
			for columnIndex, (columnSize, columnHeader) in enumerate(izip(cellSizes, headers)):
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

		for columnSpec, cellSize in izip(self.columns, cellSizes):
			columnSpec["length"] = cellSize


	def append(self, sim):
		self.nLines += 1

		if self.nLines % self.iterFreq != 0:
			return

		if self.nLines % self.headerFreq == 0:
			self.printHeaders()

		for iColumn in xrange(len(self.columns)):
			column = self.columns[iColumn]

			if iColumn > 0:
				sys.stdout.write(SPACER)

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

			sys.stdout.write(("%" + str(column["length"]) + column["format"]) % value)

		sys.stdout.write("\n")


	def finalize(self, sim):
		# Print summary
		sys.stdout.write("\n")
		sys.stdout.write("Simulation finished:\n")

		simTime = sim.time()

		# Length
		h = np.floor(simTime // 3600)
		m = np.floor((simTime % 3600) // 60)
		s = simTime % 60
		sys.stdout.write(" - Length: %d:%02d:%02.0f\n" % (h, m, s))

		# Runtime
		diff = time.time() - self.startTime
		m, s = divmod(diff, 60)
		h, m = divmod(m, 60)
		sys.stdout.write(" - Runtime: %d:%02d:%02.0f\n" % (h, m, s))

		# New line
		sys.stdout.write("\n")

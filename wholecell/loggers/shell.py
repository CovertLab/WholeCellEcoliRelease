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

import wholecell.loggers.logger

class Shell(wholecell.loggers.logger.Logger):
	""" Shell """

	def __init__(self):
		self.iterFreq = 1
		self.headerFreq = 50

	def initialize(self, sim):
		# Array of columns
		self.columns = [
			{"header": "Time (s)", "target": "Simulation", "property": "time", "length": 8, "format": "d", "sum": False},
			{"header": "Dry mass (fg)", "target": "Listener:Mass", "property": "cellDry", "length": 13, "format": ".2f", "sum": False},
			{"header": "fold", "target": "Listener:Mass", "property": "cellDryFoldChange", "length": 5, "format": ".3f", "sum": False},
			{"header": "expected", "target": "Listener:Mass", "property": "expectedFoldChange", "length": 8, "format": ".3f", "sum": False},
			{"header": "Growth (fg/s)", "target": "Listener:Mass", "property": "growth", "length": 13, "format": ".4f", "sum": False},
			{"header": "Protein frac", "target": "Listener:Mass", "property": "proteinFraction", "length": 12, "format": ".3f", "sum": False},
			{"header": "fold", "target": "Listener:Mass", "property": "proteinFoldChange", "length": 5, "format": ".3f", "sum": False},
			{"header": "RNA frac", "target": "Listener:Mass", "property": "rnaFraction", "length": 8, "format": ".3f", "sum": False},
			{"header": "fold", "target": "Listener:Mass", "property": "rnaFoldChange", "length": 5, "format": ".3f", "sum": False},

			]

		# Collect Metadata
		self.nLines = -1
		self.startTime = time.time()

		# Print initial state
		self.append(sim)


	def printHeaders(self):
		# Print headers

		if self.nLines > 0:			
			for iColumn in xrange(len(self.columns)):
				if iColumn > 0:
					sys.stdout.write("  ")

				sys.stdout.write(("%" + str(self.columns[iColumn]["length"]) + "s") % ("=" * self.columns[iColumn]["length"]))

			sys.stdout.write("\n")

		for iColumn in xrange(len(self.columns)):
			if iColumn > 0:
				sys.stdout.write("  ")

			sys.stdout.write(("%" + str(self.columns[iColumn]["length"]) + "s") % self.columns[iColumn]["header"])

		sys.stdout.write("\n")

		for iColumn in xrange(len(self.columns)):
			if iColumn > 0:
				sys.stdout.write("  ")

			sys.stdout.write(("%" + str(self.columns[iColumn]["length"]) + "s") % ("=" * self.columns[iColumn]["length"]))

		sys.stdout.write("\n")


	def append(self, sim):
		self.nLines += 1

		if self.nLines % self.iterFreq != 0:
			return

		if self.nLines % self.headerFreq == 0:
			self.printHeaders()

		for iColumn in xrange(len(self.columns)):
			column = self.columns[iColumn]

			if iColumn > 0:
				sys.stdout.write("  ")

			if column["target"] == "Simulation":
				target = sim

			else:
				targetType, targetName = column["target"].split(":")
				
				target = {
					"State":sim.states,
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

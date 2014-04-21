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

	def initialize(self, sim):
		# Array of columns
		self.columns = [
			{"header": "Time (s)", "state": "Simulation", "property": "time", "length": 8, "format": "d", "sum": False},
			{"header": "Mass (fg)", "state": "Mass", "property": "cell", "length": 9, "format": ".2f", "sum": True},
			{"header": "Growth (fg/s)", "state": "Mass", "property": "growth", "length": 13, "format": ".2f", "sum": False}
			]

		# Collect Metadata
		self.nLines = -1
		self.startTime = time.time()

		# Print headers
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

		# Print initial state
		self.append(sim)


	def append(self, sim):
		self.nLines += 1

		if self.nLines % self.iterFreq != 0:
			return

		for iColumn in xrange(len(self.columns)):
			column = self.columns[iColumn]

			if iColumn > 0:
				sys.stdout.write("  ")

			if column["state"] == "Simulation":
				target = sim

			else:
				target = sim.states[column["state"]]

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

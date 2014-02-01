#!/usr/bin/env python

"""
Shell

Prints a very brief summary of a whole-cell simulation to standard output

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

import time
import sys
import numpy

import wholecell.sim.logger.Logger

class Shell(wholecell.sim.logger.Logger.Logger):
	""" Shell """

	def __init__(self):
		self.iterFreq = 1

	def initialize(self, sim):
		# Array of columns
		self.columns = [
			{"header": "Time (s)", "state": "Time", "property": "value", "length": 8, "format": "d", "sum": False},
			{"header": "Mass (fg)", "state": "Mass", "property": "cell", "length": 9, "format": ".2f", "sum": True},
			{"header": "Growth (fg/s)", "state": "Mass", "property": "growth", "length": 13, "format": ".2f", "sum": False}
			]

		# Collect Metadata
		self.nLines = -1
		self.startTime = time.time()

		# Print headers
		for iColumn in xrange(len(self.columns)):
			self.columns[iColumn]["stateIdx"] = sim.getStateIndex(self.columns[iColumn]["state"])
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

			val = getattr(sim.states[column["stateIdx"]], column["property"])

			if column["sum"]:
				val = numpy.sum(val)

			sys.stdout.write(("%" + str(column["length"]) + column["format"]) % val)

		sys.stdout.write("\n")


	def finalize(self, sim):
		# Print summary
		sys.stdout.write("\n")
		sys.stdout.write("Simulation finished:\n")

		# Length
		h = numpy.floor(sim.getState("Time").value // 3600)
		m = numpy.floor((sim.getState("Time").value % 3600) // 60)
		s = sim.getState("Time").value % 60
		sys.stdout.write(" - Length: %d:%02d:%02.0f\n" % (h, m, s))

		# Runtime
		diff = time.time() - self.startTime
		m, s = divmod(diff, 60)
		h, m = divmod(m, 60)
		sys.stdout.write(" - Runtime: %d:%02d:%02.0f\n" % (h, m, s))

		# New line
		sys.stdout.write("\n")

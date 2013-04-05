#!/usr/bin/env python

"""
RnaMaturation

RNA maturation sub-model. Encodes molecular simulation of RNA maturation.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

import numpy

import wholecell.sim.process.Process

class RnaMaturation(wholecell.sim.process.Process.Process):
	""" RnaMaturation """

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "RnaMaturation",
			"name": "RNA maturation"
		}

		# References to states
		self.nascentRna = None
		self.matureRna = None

		super(RnaMaturation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(RnaMaturation, self).initialize(sim, kb)

		self.nascentRna = sim.getState("MoleculeCounts").addPartition(self, [x["id"] + ":nascent[c]" for x in kb.rnas], self.calcReqNascentRna)
		self.matureRna = sim.getState("MoleculeCounts").addPartition(self, [x["id"] + ":mature[c]" for x in kb.rnas], self.calcReqMatureRna)

	# Calculate needed nascent RNA
	def calcReqNascentRna(self):
		return numpy.ones(self.nascentRna.fullCounts.shape)

	# Calculate needed mature RNA
	def calcReqMatureRna(self):
		return numpy.zeros(self.matureRna.fullCounts.shape)

	# Calculate temporal evolution
	def evolveState(self):
		self.matureRna.counts += self.nascentRna.counts
		self.nascentRna.counts[:] = 0
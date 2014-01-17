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

		super(RnaMaturation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(RnaMaturation, self).initialize(sim, kb)

		mc = sim.getState('MoleculeCounts')

		self.nascentRna = mc.addPartition(self, [x["id"] + ":nascent[c]" for x in kb.rnas], self.calcReqNascentRna)
		self.matureRna = mc.addPartition(self, [x["id"] + ":mature[c]" for x in kb.rnas], self.calcReqMatureRna)


	def calcReqNascentRna(self, request):
		request.countsBulkIs(1)


	def calcReqMatureRna(self, request):
		request.countsBulkIs(0)


	# Calculate temporal evolution
	def evolveState(self):
		self.matureRna.countsBulkInc(self.nascentRna.countsBulk())
		self.nascentRna.countsBulkIs(0)

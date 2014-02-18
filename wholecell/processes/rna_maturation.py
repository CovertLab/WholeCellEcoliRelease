#!/usr/bin/env python

"""
RnaMaturation

RNA maturation sub-model. Encodes molecular simulation of RNA maturation.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

import numpy

import wholecell.processes.process

class RnaMaturation(wholecell.processes.process.Process):
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

		mc = sim.states['MoleculeCounts']

		nascentRnaIds = [x["id"] + ":nascent[c]" for x in kb.rnas]
		matureRnaIds = [x["id"] + ":mature[c]" for x in kb.rnas]

		self.mcPartition.initialize(nascentRnaIds + matureRnaIds)
		
		self.mcPartition.nascentRna = self.mcPartition.countsBulkViewNew(nascentRnaIds)
		self.mcPartition.matureRna = self.mcPartition.countsBulkViewNew(matureRnaIds)


	def requestMoleculeCounts(self):
		self.mcPartition.nascentRna.countsBulkIs(1)


	# Calculate temporal evolution
	def evolveState(self):
		self.mcPartition.matureRna.countsBulkInc(
			self.mcPartition.nascentRna.countsBulk()
			)

		self.mcPartition.nascentRna.countsBulkIs(0)

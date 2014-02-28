#!/usr/bin/env python

"""
RnaMaturation

RNA maturation sub-model. Encodes molecular simulation of RNA maturation.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

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

		mc = sim.states['BulkCounts']

		nascentRnaIds = [x["id"] + ":nascent[c]" for x in kb.rnas]
		matureRnaIds = [x["id"] + ":mature[c]" for x in kb.rnas]

		self.bulkCountsPartition.initialize(nascentRnaIds + matureRnaIds)
		
		self.bulkCountsPartition.nascentRna = self.bulkCountsPartition.countsBulkViewNew(nascentRnaIds)
		self.bulkCountsPartition.matureRna = self.bulkCountsPartition.countsBulkViewNew(matureRnaIds)


	def requestBulkCounts(self):
		self.bulkCountsPartition.nascentRna.countsBulkIs(1)


	# Calculate temporal evolution
	def evolveState(self):
		self.bulkCountsPartition.matureRna.countsBulkInc(
			self.bulkCountsPartition.nascentRna.countsBulk()
			)

		self.bulkCountsPartition.nascentRna.countsBulkIs(0)

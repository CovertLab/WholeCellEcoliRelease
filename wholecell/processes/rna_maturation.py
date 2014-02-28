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

		mc = sim.states['BulkMolecules']

		nascentRnaIds = [x["id"] + ":nascent[c]" for x in kb.rnas]
		matureRnaIds = [x["id"] + ":mature[c]" for x in kb.rnas]

		self.bulkMoleculesPartition.initialize(nascentRnaIds + matureRnaIds)
		
		self.bulkMoleculesPartition.nascentRna = self.bulkMoleculesPartition.countsBulkViewNew(nascentRnaIds)
		self.bulkMoleculesPartition.matureRna = self.bulkMoleculesPartition.countsBulkViewNew(matureRnaIds)


	def requestBulkMolecules(self):
		self.bulkMoleculesPartition.nascentRna.countsBulkIs(1)


	# Calculate temporal evolution
	def evolveState(self):
		self.bulkMoleculesPartition.matureRna.countsBulkInc(
			self.bulkMoleculesPartition.nascentRna.countsBulk()
			)

		self.bulkMoleculesPartition.nascentRna.countsBulkIs(0)

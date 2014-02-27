#!/usr/bin/env python

"""
ProteinMaturation

Protein maturation sub-model. Encodes molecular simulation of protein maturation: processing, localization

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

import wholecell.processes.process

class ProteinMaturation(wholecell.processes.process.Process):
	""" ProteinMaturation """

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "ProteinMaturation",
			"name": "ProteinMaturation",
		}

		super(ProteinMaturation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(ProteinMaturation, self).initialize(sim, kb)

		mc = sim.states["MoleculeCounts"]

		monomers = [x for x in kb.proteins if len(x["composition"]) == 0]

		nascentMonomerIds = [x["id"] + ":nascent[c]" for x in monomers]
		matureMonomerIds = [x["id"] + ":mature[" + x["location"] + "]" for x in monomers]

		self.mcPartition.initialize(nascentMonomerIds + matureMonomerIds)

		self.mcPartition.nascentMonomers = self.mcPartition.countsBulkViewNew(
			nascentMonomerIds)

		self.mcPartition.matureMonomers = self.mcPartition.countsBulkViewNew(
			matureMonomerIds)


	def requestMoleculeCounts(self):
		self.mcPartition.nascentMonomers.countsBulkIs(1)


	# Calculate temporal evolution
	def evolveState(self):
		self.mcPartition.matureMonomers.countsBulkInc(
			self.mcPartition.nascentMonomers.countsBulk())

		self.mcPartition.nascentMonomers.countsBulkIs(0)

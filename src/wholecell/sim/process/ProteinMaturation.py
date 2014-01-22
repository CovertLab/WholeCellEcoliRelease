#!/usr/bin/env python

"""
ProteinMaturation

Protein maturation sub-model. Encodes molecular simulation of protein maturation: processing, localization

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

import numpy

import wholecell.sim.process.Process

class ProteinMaturation(wholecell.sim.process.Process.Process):
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

		monomers = [x for x in kb.proteins if len(x["composition"]) == 0]

		self.nascentProteinMonomerPartition = sim.getState("MoleculeCounts").addPartition(self,
			[x["id"] + ":nascent[c]" for x in monomers],
			self.calcReqNascentProteinMonomer)

		self.matureProteinMonomerPartition = sim.getState("MoleculeCounts").addPartition(self,
			[x["id"] + ":mature[" + x["location"] + "]" for x in monomers],
			self.calcReqMatureProteinMonomer)


	# Calculate needed proteins
	def calcReqNascentProteinMonomer(self, request):
		request.countsBulkIs(1)


	# Calculate needed proteins
	def calcReqMatureProteinMonomer(self, request):
		request.countsBulkIs(0)


	# Calculate temporal evolution
	def evolveState(self):
		self.matureProteinMonomerPartition.countsBulkInc(
			self.nascentProteinMonomerPartition.countsBulk())

		self.nascentProteinMonomerPartition.countsBulkIs(0)

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

		# References to states
		self.nascentProteinMonomer = None
		self.matureProteinMonomer = None
		
		super(ProteinMaturation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(ProteinMaturation, self).initialize(sim, kb)

		monomers = [x for x in kb.proteins if len(x["composition"]) == 0]

		self.nascentProteinMonomer = sim.getState("MoleculeCounts").addPartition(self,
			[x["id"] + ":nascent[c]" for x in monomers],
			_calcReqNascentProteinMonomer)

		self.matureProteinMonomer = sim.getState("MoleculeCounts").addPartition(self,
			[x["id"] + ":mature[" + x["location"] + "]" for x in monomers],
			_calcReqMatureProteinMonomer)

	# Calculate temporal evolution
	def evolveState(self):
		self.matureProteinMonomer.countsBulkInc(
			self.nascentProteinMonomer.countsBulk())
		self.nascentProteinMonomer.countsBulkIs(0)


# Calculate needed proteins
def _calcReqNascentProteinMonomer(partition):
	return numpy.ones_like(partition.countsBulk())

# Calculate needed proteins
def _calcReqMatureProteinMonomer(partition):
	return numpy.zeros_like(partition.countsBulk())
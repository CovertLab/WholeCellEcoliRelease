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

		monomers = [x for x in kb.proteins if x["monomer"] == True]

		self.nascentProteinMonomer = sim.getState("MoleculeCounts").addPartition(self,
			[x["id"] + ":nascent[c]" for x in monomers],
			self.calcReqNascentProteinMonomer)

		self.matureProteinMonomer = sim.getState("MoleculeCounts").addPartition(self,
			[x["id"] + ":mature[" + x["location"] + "]" for x in monomers],
			self.calcReqMatureProteinMonomer)

	# Calculate needed proteins
	def calcReqNascentProteinMonomer(self):
		return numpy.ones(self.nascentProteinMonomer.fullCounts.shape)

	# Calculate needed proteins
	def calcReqMatureProteinMonomer(self):
		return numpy.zeros(self.matureProteinMonomer.fullCounts.shape)

	# Calculate temporal evolution
	def evolveState(self):
		self.matureProteinMonomer.counts += self.nascentProteinMonomer.counts
		self.nascentProteinMonomer.counts[:] = 0
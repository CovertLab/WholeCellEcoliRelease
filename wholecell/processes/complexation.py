#!/usr/bin/env python

"""
Complexation

Macromolecular complexation sub-model. Encodes molecular simulation of macromolecular complexation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

import numpy as np

import wholecell.processes.process

class Complexation(wholecell.processes.process.Process):
	""" Complexation """

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "Complexation",
			"name": "Macromolecular complexation",
		}

		# References to states
		self.subunit = None
		self.complex = None
		
		super(Complexation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Complexation, self).initialize(sim, kb)

		# # Complex
		# complexes = [x for x in kb.proteins if x["monomer"] == False and x["formationProcess"] == self.meta["id"]]
		# self.complex = sim.getState("MoleculeCounts").addPartition(self,
		# 	[x["id"] + ":mature[" + x["location"] + "]" for x in complexes],
		# 	self.calcReqComplex)

		# # Subunits
		# subunits = []
		# for c in complexes:
		# 	subunits.extend([x for x in c["composition"] if x["coeff"] < 0])
		# subIdComps = list(set([x["molecule"] + ":mature[" + x["location"] + "]" for x in subunits]))

		# self.subunit = sim.getState("MoleculeCounts").addPartition(self, subIdComps, self.calcReqSubunit)

		# tmpSMat = []
		# for iComplex in xrange(len(complexes)):
		# 	c = complexes[iComplex]
		# 	for s in c["composition"]:
		# 		if s["coeff"] > 0:
		# 			continue
		# 		tmpSMat.append([s["molecule"] + ":mature[" + s["location"] + "]", iComplex, -s["coeff"]])

		# subIdx = self.subunit.getIndex([x[0] for x in tmpSMat])[0]
		# self.sMat = np.zeros((len(subIdComps), len(complexes)))
		# self.sMat[subIdx, np.array([x[1] for x in tmpSMat])] = np.array([x[2] for x in tmpSMat])

	# Calculate needed proteins (subunits)
	def calcReqSubunit(self):
		return np.ones(self.subunit.fullCounts.shape)

	# Calculate needed proteins (complexes)
	def calcReqComplex(self):
		return np.zeros(self.complex.fullCounts.shape)

	# Calculate temporal evolution
	def evolveState(self):
		return
		self.subunit.counts, self.complex.counts = self.calcNewComplexes(self.subunit.counts, self.complex.counts, 1)

	# Gillespie-like algorithm
	def calcNewComplexes(self, subunits, complexes, leap):
		return
		import warnings
		warnings.simplefilter("ignore", RuntimeWarning)	# Supress warnings about divide by zero
		# TODO: Implement tau leaping correctly
		while True:
			# Calculate rates
			rates = np.floor(np.nanmin(			# For each complex, set its rate to be that of its limiting subunit
				subunits.reshape((-1,1)) / self.sMat	# Note: np's broadcasting mechanisms tile subunits to be the shape of sMat
				, axis = 0))
			totRate = np.sum(rates)
			rates /= totRate
			# i += 1
			# if i > 50:
			# 	import ipdb
			# 	ipdb.set_trace()
			if totRate <= 0:
				break

			# Check if sufficient metabolic resources to make protein
			newCnts = self.randStream.mnrnd(leap, rates)
			if np.any(np.dot(self.sMat, newCnts) > subunits):
				break

			# Update subunits
			subunits -= np.dot(self.sMat, newCnts)

			# Increment complexes
			complexes += newCnts
		return subunits, complexes

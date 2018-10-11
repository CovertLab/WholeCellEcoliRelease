#!/usr/bin/env python

"""
Complexation

Macromolecular complexation sub-model. Encodes molecular simulation of macromolecular complexation

TODO:
- allow for shuffling when appropriate (maybe in another process)
- handle protein complex dissociation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices

# Maximum unsigned int value + 1 for randint() to seed srand from C stdlib
RAND_MAX = 2**32

class Complexation(wholecell.processes.process.Process):
	""" Complexation """

	_name = "Complexation"

	# Constructor
	def __init__(self):

		super(Complexation, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(Complexation, self).initialize(sim, sim_data)

		# Create matrices and vectors that describe reaction stoichiometries 
		self.stoichMatrix = sim_data.process.complexation.stoichMatrix().astype(np.int64, order = "F")

		self.prebuiltMatrices = mccBuildMatrices(self.stoichMatrix)

		self.product_indices = [idx for idx in np.where(np.any(self.stoichMatrix > 0, axis=1))[0]]

		# Build views
		moleculeNames = sim_data.process.complexation.moleculeNames
		self.molecules = self.bulkMoleculesView(moleculeNames)


	def calculateRequest(self):
		moleculeCounts = self.molecules.total()

		# Macromolecule complexes are requested
		updatedMoleculeCounts, complexationEvents = mccFormComplexesWithPrebuiltMatrices(
			moleculeCounts,
			self.randomState.randint(RAND_MAX),
			self.stoichMatrix,
			*self.prebuiltMatrices
			)

		self.molecules.requestIs(np.fmax(moleculeCounts - updatedMoleculeCounts, 0))


	def evolveState(self):
		moleculeCounts = self.molecules.counts()

		# Macromolecule complexes are formed from their subunits
		updatedMoleculeCounts, complexationEvents = mccFormComplexesWithPrebuiltMatrices(
			moleculeCounts,
			self.randomState.randint(RAND_MAX),
			self.stoichMatrix,
			*self.prebuiltMatrices
			)

		self.molecules.countsIs(updatedMoleculeCounts)

		# Write outputs to listeners
		self.writeToListener("ComplexationListener", "complexationEvents", complexationEvents)

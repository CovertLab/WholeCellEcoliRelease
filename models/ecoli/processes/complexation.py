#!/usr/bin/env python

"""
Complexation

Macromolecular complexation sub-model. Encodes molecular simulation of macromolecular complexation

TODO:
- allow for shuffling when appropriate (maybe in another process)
- handle complex degradation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.mc_complexation import mccBuildMatrices, mccFormComplexesWithPrebuiltMatrices

class Complexation(wholecell.processes.process.Process):
	""" Complexation """

	_name = "Complexation"

	# Constructor
	def __init__(self):

		super(Complexation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(Complexation, self).initialize(sim, kb)

		# Create matrices and vectors

		self.stoichMatrix = kb.process.complexation.stoichMatrix().astype(np.int64, order = "F")

		self.prebuiltMatrices = mccBuildMatrices(self.stoichMatrix)

		# Build views

		moleculeNames = kb.process.complexation.moleculeNames

		self.molecules = self.bulkMoleculesView(moleculeNames)


	def calculateRequest(self):
		moleculeCounts = self.molecules.total()

		updatedMoleculeCounts = mccFormComplexesWithPrebuiltMatrices(
			moleculeCounts,
			self.seed,
			self.stoichMatrix,
			*self.prebuiltMatrices
			)

		self.molecules.requestIs(np.fmax(moleculeCounts - updatedMoleculeCounts, 0))


	def evolveState(self):
		moleculeCounts = self.molecules.counts()

		updatedMoleculeCounts = mccFormComplexesWithPrebuiltMatrices(
			moleculeCounts,
			self.seed,
			self.stoichMatrix,
			*self.prebuiltMatrices
			)

		self.molecules.countsIs(updatedMoleculeCounts)

#!/usr/bin/env python

"""
Complexation

Macromolecular complexation sub-model. Encodes molecular simulation of macromolecular complexation

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.mc_complexation import monteCarloComplexation

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

		self.stoichMatrix = kb.complexationStoichMatrix().astype(np.int64, order = "F")

		# Build views

		moleculeNames = kb.complexationMoleculeNames
		subunitNames = kb.complexationSubunitNames
		# complexNames = kb.complexationComplexNames

		self.molecules = self.bulkMoleculesView(moleculeNames)
		self.subunits = self.bulkMoleculesView(subunitNames)

		self.moleculeNames = moleculeNames


	def calculateRequest(self):
		self.subunits.requestAll()


	def evolveState(self):
		moleculeCounts = self.molecules.counts()

		updatedMoleculeCounts = monteCarloComplexation(
			moleculeCounts,
			self.stoichMatrix,
			self.seed
			)

		self.molecules.countsIs(updatedMoleculeCounts)

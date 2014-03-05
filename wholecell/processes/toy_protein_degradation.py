#!/usr/bin/env python

"""
ToyProteinDegradation

Degrages imaginary RNA polymerases only if they are unbound.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/28/2014
"""

import numpy as np

import wholecell.processes.process

class ToyProteinDegradation(wholecell.processes.process.Process):
	""" ToyProteinDegradation """

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "ToyProteinDegradation",
			"name": "ToyProteinDegradation",
		}

		self.degradationProbability = 0.1

		super(ToyProteinDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(ToyProteinDegradation, self).initialize(sim, kb)

		uniqueMolecules = sim.states["UniqueMolecules"]
		self.container = uniqueMolecules._container

	def requestMoleculeCounts(self):
		# TODO: Needs to request bulk RNA polymerase molecules
		pass

	# Calculate temporal evolution
	def evolveState(self):
		unbound = self.container.evaluateQuery('RNA polymerase', boundToChromosome = ('==', False))
		toDelete = []
		for molecule in unbound:
			if self.randStream.rand <= self.degradationProbability:
				toDelete.append(molecule)

		for molecule in toDelete:
			self.container.moleculeDel(molecule)

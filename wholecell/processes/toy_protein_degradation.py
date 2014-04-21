#!/usr/bin/env python

"""
ToyProteinDegradation

Degrages imaginary RNA polymerases only if they are unbound.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/28/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class ToyProteinDegradation(wholecell.processes.process.Process):
	""" ToyProteinDegradation """

	_name = "ToyProteinDegradation"

	# Constructor
	def __init__(self):
		self.degradationProbability = 0.1

		super(ToyProteinDegradation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(ToyProteinDegradation, self).initialize(sim, kb)

		self.unboundRNApoly = self.uniqueMoleculesView('RNA polymerase',
			boundToChromosome = ('==', False))


	def calculateRequest(self):
		nDegrade = self.unboundRNApoly.total() * self.degradationProbability

		self.unboundRNApoly.requestIs(nDegrade)


	# Calculate temporal evolution
	def evolveState(self):
		unbound = self.unboundRNApoly.molecules()

		self.unboundRNApoly.moleculesDel(unbound)

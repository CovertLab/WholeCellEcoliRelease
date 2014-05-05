#!/usr/bin/env python

"""
MetabolismFba

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/5/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

UNCONSTRAINED_FLUX_VALUE = 1000.0

class MetabolismFba(wholecell.processes.process.Process):
	""" MetabolismFba """

	_name = "MetabolismFba"

	# Constructor
	def __init__(self):
		self.lowerBounds = None
		self.upperBounds = None
		self.stoichMatrix = None
		self.objective = None

		super(MetabolismFba, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(MetabolismFba, self).initialize(sim, kb)

		# TODO: see if biomass reaction is in some way already in the table
		
		wildtypeIds = kb.wildtypeBiomass['metaboliteId']
		wildtypeBiomassReaction = kb.wildtypeBiomass['biomassFlux'].magnitude

		# Must add one entry for the biomass reaction

		stoichMatrix = kb.metabolismStoichMatrix

		self.stoichMatrix = np.hstack([
			stoichMatrix,
			np.zeros((stoichMatrix.shape[0], 1))
			])

		self.lowerBounds = np.zeros(self.stoichMatrix.shape[1])
		self.lowerBounds[kb.metabolismReversibleReactions] = -UNCONSTRAINED_FLUX_VALUE

		self.upperBounds = UNCONSTRAINED_FLUX_VALUE * np.ones(self.stoichMatrix.shape[1])

		indexes = [kb.metabolismMoleculeNames.index(moleculeName) for moleculeName in wildtypeIds]

		self.stoichMatrix[indexes, -1] = wildtypeBiomassReaction

		self.objective = np.zeros(self.stoichMatrix.shape[1])
		self.objective[-1] = -1 # TODO: check signs

		self.molecules = self.bulkMoleculesView(kb.metabolismMoleculeNames)


	def calculateRequest(self):
		import ipdb; ipdb.set_trace()

		# TODO: how to generate request?


	# Calculate temporal evolution
	def evolveState(self):
		print 'evolve'

#!/usr/bin/env python

"""
MetabolismFba

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/5/2014
"""

from __future__ import division

import warnings

import numpy as np

import wholecell.processes.process

UNCONSTRAINED_FLUX_VALUE = 10000.0

class MetabolismFba(wholecell.processes.process.Process):
	""" MetabolismFba """

	_name = "MetabolismFba"

	# Construct object graph
	def initialize(self, sim, kb):
		super(MetabolismFba, self).initialize(sim, kb)

		# TODO: see if biomass reaction is in some way already in the table
		
		wildtypeIds = kb.wildtypeBiomass['metaboliteId']
		self.biomassReaction = kb.wildtypeBiomass['biomassFlux'].magnitude * 1e-3 * kb.nAvogadro.to('1 / mole').magnitude * kb.avgCellDryMassInit.to('g').magnitude

		# Must add one entry for the biomass reaction

		self.stoichMatrix = np.hstack([
			kb.metabolismStoichMatrix,
			np.zeros((kb.metabolismStoichMatrix.shape[0], 1))
			])

		self.nFluxes = self.stoichMatrix.shape[1]

		self.reversibleReactions = kb.metabolismReversibleReactions

		indexes = [kb.metabolismMoleculeNames.index(moleculeName) for moleculeName in wildtypeIds]

		self.stoichMatrix[indexes, -1] = -self.biomassReaction

		self.objective = np.zeros(self.nFluxes)
		self.objective[-1] = -1 # TODO: check signs

		self.exchangeReactionIndexes = kb.metabolismExchangeReactionIndexes

		self.exchangeReactionMolecules = kb.metabolismExchangeReactionMolecules

		self.exchangedMolecules = self.bulkMoleculesView(
			self.exchangeReactionMolecules
			)

		self.biomassMolecules = self.bulkMoleculesView(wildtypeIds)


	def calculateRequest(self):
		# TODO: eliminate the LP solving here
		self.exchangedMolecules.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# Set up LP

		lowerBounds = np.zeros(self.nFluxes)
		lowerBounds[self.reversibleReactions] = -UNCONSTRAINED_FLUX_VALUE

		upperBounds = np.empty(self.nFluxes)
		upperBounds.fill(UNCONSTRAINED_FLUX_VALUE)

		# for moleculeIndex in np.where(self.exchangedMolecules.counts() == 0)[0]:
		# 	warnings.warn('Unavailable exchange reaction molecule: {}'.format(
		# 		self.exchangeReactionMolecules[moleculeIndex]
		# 		))

		upperBounds[self.exchangeReactionIndexes] = self.exchangedMolecules.counts() / self.timeStepSec
		
		# upperBounds[self.exchangeReactionIndexes] = np.fmax(
		# 	UNCONSTRAINED_FLUX_VALUE,
		# 	self.exchangedMolecules.counts() / self.timeStepSec
		# 	)

		fluxes, status = fba(self.stoichMatrix, lowerBounds, upperBounds, self.objective)

		if status != "optimal":
			warnings.warn("Linear programming did not converge")

		if np.any(np.abs(fluxes) == UNCONSTRAINED_FLUX_VALUE):
			warnings.warn("Reaction fluxes reached 'unconstrained' boundary")

		exchangeUsage = -(fluxes[self.exchangeReactionIndexes] * self.timeStepSec).astype(np.int)

		biomassProduction = (fluxes[-1] * self.timeStepSec * self.biomassReaction).astype(np.int)

		# import ipdb; ipdb.set_trace()

		self.exchangedMolecules.countsDec(exchangeUsage)

		self.biomassMolecules.countsInc(biomassProduction)

import cvxopt.solvers
from cvxopt import matrix, sparse, spmatrix

def fba(stoichiometricMatrix, lowerBounds, upperBounds, objective):
	cvxopt.solvers.options["LPX_K_MSGLEV"] = 0

	nNodes, nEdges = stoichiometricMatrix.shape

	A = sparse(matrix(stoichiometricMatrix)) # NOTE: I don't know if this actually helps the solver
	h = matrix(np.concatenate([upperBounds, -lowerBounds], axis = 0))
	f = matrix(objective)

	b = matrix(np.zeros(nNodes))
	G = spmatrix(
		[1]*nEdges + [-1]*nEdges,
		np.arange(2*nEdges),
		np.tile(np.arange(nEdges), 2)
		)

	solution = cvxopt.solvers.lp(f, G, h, A = A, b = b, solver = 'glpk')

	fluxes = np.array(solution['x']).flatten()

	status = solution["status"]

	return fluxes, status

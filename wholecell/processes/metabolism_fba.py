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

REQUEST_SIMPLE = False # if True, doesn't run FBA twice, but requests all possible metabolites

# TODO: better requests
# TODO: flexFBA etc
# TODO: explore dynamic biomass objectives
# TODO: dark energy accounting
# TODO: eliminate futile cycles
# TODO: cache FBA vectors/matrices instead of rebuilding
# TODO: media exchange constraints

class MetabolismFba(wholecell.processes.process.Process):
	""" MetabolismFba """

	_name = "MetabolismFba"

	# Construct object graph
	def initialize(self, sim, kb):
		super(MetabolismFba, self).initialize(sim, kb)
		
		self.biomassIds = kb.wildtypeBiomass['metaboliteId']
		self.biomassReaction = ( # TODO: validate this math
			kb.wildtypeBiomass['biomassFlux'].magnitude
			* 1e-3
			* kb.nAvogadro.to('1 / mole').magnitude
			* kb.avgCellDryMassInit.to('g').magnitude
			)

		# Must add one entry for the biomass reaction

		self.stoichMatrix = np.hstack([
			kb.metabolismStoichMatrix,
			np.zeros((kb.metabolismStoichMatrix.shape[0], 1))
			])

		self.nFluxes = self.stoichMatrix.shape[1]

		self.reversibleReactions = np.hstack([
			kb.metabolismReversibleReactions,
			np.array([False], np.bool)
			])

		indexes = [kb.metabolismMoleculeNames.index(moleculeName) for moleculeName in self.biomassIds]

		self.stoichMatrix[indexes, -1] = -self.biomassReaction

		self.objective = np.zeros(self.nFluxes)
		self.objective[-1] = -1

		self.mediaExchangeMoleculeNames = kb.metabolismMediaExchangeReactionNames
		self.mediaExchangeIndexes = kb.metabolismMediaExchangeReactionIndexes

		self.sinkExchangeMoleculeNames = kb.metabolismSinkExchangeReactionNames
		self.sinkExchangeIndexes = kb.metabolismSinkExchangeReactionIndexes

		self.internalExchangeMoleculeNames = kb.metabolismInternalExchangeReactionNames
		self.internalExchangeIndexes = kb.metabolismInternalExchangeReactionIndexes

		self.biomassInternalExchangeIndexes = np.array([
			self.internalExchangeIndexes[self.internalExchangeMoleculeNames.index(moleculeId)]
			for moleculeId in self.biomassIds
			if moleculeId in self.internalExchangeMoleculeNames
			])

		# Create views

		self.biomassMolecules = self.bulkMoleculesView(self.biomassIds)

		self.sinkMolecules = self.bulkMoleculesView(self.sinkExchangeMoleculeNames)

		self.internalExchangeMolecules = self.bulkMoleculesView(self.internalExchangeMoleculeNames)


	def calculateRequest(self):
		if REQUEST_SIMPLE:
			# Request every internal metabolite
			self.internalExchangeMolecules.requestAll()

		else:
			# Compute request assuming all metabolites are allocated
			totalCounts = self.internalExchangeMolecules.total()

			fluxes = self.computeFluxes(totalCounts)

			internalUsage = np.fmax(
				(fluxes[self.internalExchangeIndexes] * self.timeStepSec).astype(np.int),
				0
				)

			self.internalExchangeMolecules.requestIs(internalUsage)

			print ', '.join(
				self.internalExchangeMoleculeNames[index]
				for index in np.where(internalUsage)[0]
				)


	# Calculate temporal evolution
	def evolveState(self):
		# Update metabolite counts based on computed fluxes

		fluxes = self.computeFluxes(self.internalExchangeMolecules.counts())

		internalUsage = (fluxes[self.internalExchangeIndexes] * self.timeStepSec).astype(np.int)
		sinkProduction = (fluxes[self.sinkExchangeIndexes] * self.timeStepSec).astype(np.int)
		biomassProduction = (fluxes[-1] * self.timeStepSec * self.biomassReaction).astype(np.int)

		self.internalExchangeMolecules.countsDec(internalUsage)
		self.sinkMolecules.countsInc(sinkProduction)
		self.biomassMolecules.countsInc(biomassProduction)


	def computeFluxes(self, internalMoleculeCounts):
		# Set up LP

		lowerBounds = np.zeros(self.nFluxes)
		lowerBounds[self.reversibleReactions] = -UNCONSTRAINED_FLUX_VALUE

		upperBounds = np.empty(self.nFluxes)
		upperBounds.fill(UNCONSTRAINED_FLUX_VALUE)

		# TODO: find actual media exchange limits
		upperBounds[self.mediaExchangeIndexes] = UNCONSTRAINED_FLUX_VALUE

		upperBounds[self.internalExchangeIndexes] = internalMoleculeCounts / self.timeStepSec

		# Forbid exchange of biomass components
		lowerBounds[self.biomassInternalExchangeIndexes] = 0
		upperBounds[self.biomassInternalExchangeIndexes] = 0

		fluxes, status = fba(self.stoichMatrix, lowerBounds, upperBounds, self.objective)

		if status != "optimal":
			warnings.warn("Linear programming did not converge")

		# if np.any(np.abs(fluxes) == UNCONSTRAINED_FLUX_VALUE):
		# 	warnings.warn("Reaction fluxes reached 'unconstrained' boundary")

		return fluxes


import cvxopt.solvers
from cvxopt import matrix, sparse, spmatrix

def fba(stoichiometricMatrix, lowerBounds, upperBounds, objective):
	# Solve the linear program:
	# 0 = Sv
	# lb <= v <= ub
	# max {f^T v}

	# Supress output
	cvxopt.solvers.options["LPX_K_MSGLEV"] = 0

	nNodes, nEdges = stoichiometricMatrix.shape

	# Create cvxopt types
	A = sparse(matrix(stoichiometricMatrix)) # NOTE: I don't know if this actually helps the solver
	h = matrix(np.concatenate([upperBounds, -lowerBounds], axis = 0))
	f = matrix(objective)

	b = matrix(np.zeros(nNodes))
	G = spmatrix(
		[1]*nEdges + [-1]*nEdges,
		np.arange(2*nEdges),
		np.tile(np.arange(nEdges), 2)
		)

	# Solve LP
	solution = cvxopt.solvers.lp(f, G, h, A = A, b = b, solver = 'glpk')

	# Parse solution
	fluxes = np.array(solution['x']).flatten()
	status = solution["status"]

	return fluxes, status

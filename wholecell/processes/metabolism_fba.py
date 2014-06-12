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

	def __init__(self):
		super(MetabolismFba, self).__init__()


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

		basicStoichMatric = kb.metabolismStoichMatrix()

		self.stoichMatrix = np.hstack([
			basicStoichMatric,
			np.zeros((basicStoichMatric.shape[0], 1))
			])

		self.nFluxes = self.stoichMatrix.shape[1]

		self.reactionIsReversible = np.append(kb.metabolismReactionIsReversible, 0)

		# TODO: bug Nick about adding -> selenocysteine to the metabolic network

		indexes = [np.where(kb.metabolismMoleculeNames == moleculeName)[0][0] for moleculeName in self.biomassIds]

		self.stoichMatrix[indexes, -1] = -self.biomassReaction

		self.objective = np.zeros(self.nFluxes)
		self.objective[-1] = -1

		self.reactionIsMediaExchange = np.append(kb.metabolismReactionIsMediaExchange, 0)
		self.reactionIsSink = np.append(kb.metabolismReactionIsSink, 0)

		# Create views

		self.molecules = self.bulkMoleculesView(kb.metabolismMoleculeNames)

		# Permanent references to evolveState variables for listener

		self.fluxes = np.zeros(self.stoichMatrix.shape[1])


	def calculateRequest(self):
		pass


	# Calculate temporal evolution
	def evolveState(self):
		# Update metabolite counts based on computed fluxes

		self.fluxes = self._computeFluxes(self.internalExchangeMolecules.counts())

		deltaMolecules = np.dot(self.stoichMatrix[:, :-1], self.fluxes[:-1])

		import ipdb; ipdb.set_trace()

		self.molecules.countsInc(deltaMolecules)


	def _computeFluxes(self):
		# TODO: constrain reactions by enzyme count

		# Set up LP

		lowerBounds = np.zeros(self.nFluxes)
		lowerBounds[self.reactionIsReversible] = -UNCONSTRAINED_FLUX_VALUE

		upperBounds = np.empty(self.nFluxes)
		upperBounds.fill(UNCONSTRAINED_FLUX_VALUE)

		# TODO: find actual media exchange limits
		upperBounds[self.mediaExchangeIndexes] = UNCONSTRAINED_FLUX_VALUE

		fluxes, status = _fba(self.stoichMatrix, lowerBounds, upperBounds, self.objective)

		if status != "optimal":
			warnings.warn("Linear programming did not converge")

		# if np.any(np.abs(fluxes) == UNCONSTRAINED_FLUX_VALUE):
		# 	warnings.warn("Reaction fluxes reached 'unconstrained' boundary")

		return fluxes


import cvxopt.solvers
from cvxopt import matrix, sparse, spmatrix

def _fba(stoichiometricMatrix, lowerBounds, upperBounds, objective):
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

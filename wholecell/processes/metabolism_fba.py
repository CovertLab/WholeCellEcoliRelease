#!/usr/bin/env python

"""
MetabolismFba

TODO:
- eliminate this process once flexFBA is online

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
			kb.wildtypeBiomass['biomassFlux'].to("mole/DCW_gram").magnitude
			* kb.nAvogadro.to('1 / mole').magnitude
			* kb.avgCellDryMassInit.to('g').magnitude
			) * np.exp(np.log(2)/3600) * (np.exp(np.log(2)/3600)-1)

		# Must add one entry for the biomass reaction

		basicStoichMatric = kb.metabolismStoichMatrix()

		self.stoichMatrix = np.hstack([
			basicStoichMatric,
			np.zeros((basicStoichMatric.shape[0], 1))
			])

		self.nFluxes = self.stoichMatrix.shape[1]

		self.reactionIsReversible = np.append(kb.metabolismReactionIsReversible, False)

		indexes = [np.where(kb.metabolismMoleculeNames == moleculeName)[0][0] for moleculeName in self.biomassIds]

		self.stoichMatrix[indexes, -1] = -self.biomassReaction

		self.objective = np.zeros(self.nFluxes)
		self.objective[-1] = -1

		self.reactionIsMediaExchange = np.append(kb.metabolismReactionIsMediaExchange, False)
		self.reactionIsSink = np.append(kb.metabolismReactionIsSink, False)

		# Create views

		self.molecules = self.bulkMoleculesView(kb.metabolismMoleculeNames)

		# Temporary attributes for debugging

		self._moleculeNames = kb.metabolismMoleculeNames


	def calculateRequest(self):
		pass


	# Calculate temporal evolution
	def evolveState(self):
		# Update metabolite counts based on computed fluxes

		fluxes = self._computeFluxes()

		deltaMolecules = np.dot(self.stoichMatrix[:, :-1], fluxes[:-1]).astype(np.int64)

		# if fluxes[-1] != 0:
		# 	print "nonzero result"

		# for moleculeIndex in np.where(deltaMolecules)[0]:
		# 	print self._moleculeNames[moleculeIndex], deltaMolecules[moleculeIndex]

		self.molecules.countsInc(deltaMolecules)


	def _computeFluxes(self):
		# TODO: constrain reactions by enzyme count

		# Set up LP

		lowerBounds = np.zeros(self.nFluxes)
		lowerBounds[self.reactionIsReversible] = -UNCONSTRAINED_FLUX_VALUE

		upperBounds = np.empty(self.nFluxes)
		upperBounds.fill(UNCONSTRAINED_FLUX_VALUE)

		# TODO: find actual media exchange limits
		# upperBounds[self.reactionIsMediaExchange] = UNCONSTRAINED_FLUX_VALUE

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

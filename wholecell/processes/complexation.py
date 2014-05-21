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

from cvxopt import glpk
from cvxopt import matrix, sparse

EVALUATE_TO_COMPLETION = False # If true, will form as many complexes as possible

# TODO: evaluate mass balance
# TODO: speed up by using only the needed portion of the stoich matrix

class Complexation(wholecell.processes.process.Process):
	""" Complexation """

	_name = "Complexation"

	# Constructor
	def __init__(self):
		
		super(Complexation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(Complexation, self).initialize(sim, kb)

		self.nSubunits = kb.complexationMatrix.shape[0]
		self.nComplexes = kb.complexationMatrix.shape[1]

		self.nMolecules = self.nSubunits + self.nComplexes
		self.nFluxes = self.nSubunits + 2 * self.nComplexes

		self.complexationStoichMatrix = np.zeros((self.nMolecules, self.nFluxes), np.float)

		self.complexationStoichMatrix[:self.nSubunits, :self.nSubunits] = np.identity(self.nSubunits)
		self.complexationStoichMatrix[:self.nSubunits, self.nSubunits:self.nSubunits + self.nComplexes] = -kb.complexationMatrix
		self.complexationStoichMatrix[self.nSubunits:self.nSubunits + self.nComplexes, self.nSubunits:self.nSubunits + self.nComplexes] = np.identity(self.nComplexes)
		self.complexationStoichMatrix[-self.nComplexes:, -self.nComplexes:] = -np.identity(self.nComplexes)

		self.baseObjective = np.hstack([np.zeros(self.nSubunits + self.nComplexes), -np.ones(self.nComplexes)])

		self.ilpb = matrix(np.zeros(self.nMolecules))

		self.ilpG = sparse(matrix(np.vstack([
			np.identity(self.nFluxes),
			-np.identity(self.nFluxes)
			])))

		self.ilpI = set(range(self.nMolecules))

		self.subunits = self.bulkMoleculesView(kb.complexationMatrixSubunitIds)

		self.complexes = self.bulkMoleculesView(kb.complexationMatrixComplexIds)


	def calculateRequest(self):
		self.subunits.requestAll()


	def evolveState(self):
		while True:
			objective = matrix(self.baseObjective) # TODO: add noise

			lb = np.zeros(self.nFluxes)

			subunitCounts = self.subunits.counts()

			maxPossibleComplexes = subunitCounts.max()

			ub = maxPossibleComplexes * np.ones(self.nFluxes)

			ub[:self.nSubunits] = subunitCounts

			h = matrix(np.hstack([ub, -lb]))

			A = sparse(matrix(self.complexationStoichMatrix))

			# Supress output
			glpk.options["LPX_K_MSGLEV"] = 0

			solution = glpk.ilp(
				objective,
				self.ilpG, h,
				A, self.ilpb,
				self.ilpI
				)

			assert solution is not None, "Failed to find an optimal solution"

			fluxes = np.array(solution[1]).flatten()

			subunitsUsed = fluxes[:self.nSubunits].astype(np.int)
			complexesMade = fluxes[-self.nComplexes:].astype(np.int)

			self.subunits.countsDec(subunitsUsed)
			self.complexes.countsInc(complexesMade)

			if not EVALUATE_TO_COMPLETION or not fluxes.any():
				break

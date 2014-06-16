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
# TODO: spedd up by performing basic heuristics (only optimizing where the 
# problem isn't degenerate)

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

		self.stoichMatrix = kb.complexationStoichMatrix()# .astype(np.int64)
		self.compositionMatrix = np.ma.masked_array(
			(-self.stoichMatrix) * (self.stoichMatrix < 0),
			(self.stoichMatrix >= 0)
			)
		# self.productMatrix = (self.stoichMatrix ) * (self.stoichMatrix  > 0)

		# Find all reactions that are unique; i.e., orthogonal to all other
		# composition vectors
		self.reactionUsesMonomersUniquely = ~(
			np.dot(self.compositionMatrix.T, self.compositionMatrix)
			* ~np.identity(self.compositionMatrix.shape[1], np.bool)
			).any(0).data # here axis 0 or 1 is valid since the matrix is symmetric

		# Build views

		moleculeNames = kb.complexationMoleculeNames
		subunitNames = kb.complexationSubunitNames
		# complexNames = kb.complexationComplexNames

		self.molecules = self.bulkMoleculesView(moleculeNames)
		self.subunits = self.bulkMoleculesView(subunitNames)


	def calculateRequest(self):
		self.subunits.requestAll()


	def evolveState(self):
		moleculeCounts = self.molecules.counts().astype(np.float64) # works with floats to use BLAS

		# Compute the max number of reactions

		# olderr = np.seterr(divide = "ignore", invalid = "ignore")
		maxReactions = (moleculeCounts // self.compositionMatrix.T).min(1).data
		# np.seterr(olderr)

		# Perform the trivial complexation reactions

		trivialReactionCounts = maxReactions * self.reactionUsesMonomersUniquely

		moleculeCounts += np.dot(self.stoichMatrix, trivialReactionCounts)

		# For the nontrivial reactions, randomly form complexes to completion

		activeReactions = ~self.reactionUsesMonomersUniquely

		while activeReactions.any():
			reactionIndex = self.randomState.choice(np.where(activeReactions)[0])

			reactionStoich = self.stoichMatrix[:, reactionIndex]

			if (-reactionStoich <= moleculeCounts).all():
				moleculeCounts += reactionStoich
				# print "performed {}".format(reactionIndex)

			else:
				activeReactions[reactionIndex] = False
				# print "culled {}".format(reactionIndex)

		# formComplexes(moleculeCounts, self.stoichMatrix, self.randomState)

		self.molecules.countsIs(moleculeCounts)


MAX_ITERATIONS = 10**9
def formComplexes(moleculeCounts, stoichiometry, randomState): # NOTE: will need to pass a seed and use a library
	moleculeCountsOut = moleculeCounts.copy()

	nReactions = stoichiometry.shape[1]

	reactionIsActive = np.ones(nReactions, np.bool) # NOTE: may need to be int-typed
	reactionCumSum = reactionIsActive.cumsum() # NOTE: may need to be explicit

	for i in xrange(MAX_ITERATIONS):
		value = randomState.rand() * reactionCumSum[-1]

		for reactionIndex in xrange(nReactions):
			if value <= reactionCumSum[reactionIndex]:
				break

		else: # NOTE: for-else probably not supported by cython
			raise Exception("Should have broken out of this for-loop")

		reactionStoich = stoichiometry[:, reactionIndex] # NOTE: matrix ordering is bad

		if (-reactionStoich <= moleculeCountsOut).all(): # NOTE: may need to be explicit
			moleculeCountsOut += reactionStoich # NOTE: may need to be explicit

		else:
			reactionIsActive[reactionIndex] = False

			if not reactionIsActive.any(): # NOTE: may need to be explicit
				break

			reactionCumSum = reactionIsActive.cumsum() # NOTE: may need to be explicit

	return moleculeCountsOut


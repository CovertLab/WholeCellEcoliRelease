
from __future__ import division

import numpy as np

from cvxopt.glpk import ilp, lp
from cvxopt import matrix, sparse, spmatrix

# TODO: line-profile
def polymerizePooledMonomers(monomerCounts, monomerDeficits, maxElongation,
		randStream, useIntegerLinearProgramming = True):
	"""polymerizePooledMonomers

	Polymerizes (i.e. for transcription, translation) using a simplified model 
	of polymer structure; the total counts of each monomer must be satisfied, 
	but the order of addition to the nascent polymer does not matter.

	Arguments:

	monomerCounts, a vector of counts of each monomer
	monomerDeficits, a matrix of monomers required (nPolymers) x (nMonomers)
	maxElongation, a scalar value that defines the maximum elongation rate
	randStream, a randStream object
	useIntegerLinearProgramming, a flag that indicates that the solver should
		use ILP over LP

	ILP may produce better results but runs slower.
	"""

	# (I)LP algorithm:
	# Choose v such that
	# b = Av
	# v < Gh
	# min f^T v
	# (all v integer-valued)

	# Build the A matrix

	nMonomers = monomerDeficits.shape[1]
	nPolymers = monomerDeficits.shape[0]

	nNodes = nMonomers + 2*nPolymers # i.e. rows
	nEdges = nMonomers * (1 + nPolymers) + nPolymers + 1 + nPolymers # i.e. columns

	nNonzeroEntries = nMonomers + 2*nMonomers*nPolymers + 5*nPolymers

	values = np.empty(nNonzeroEntries)
	rowIndexes = np.empty(nNonzeroEntries, np.int)
	colIndexes = np.empty(nNonzeroEntries, np.int)

	values[:nMonomers] = 1
	values[nMonomers:nMonomers*(1+nPolymers)] = -1
	values[nMonomers*(1+nPolymers):nMonomers*(1+2*nPolymers)] = 1
	values[nMonomers*(1+2*nPolymers):nMonomers*(1+2*nPolymers)+nPolymers] = -1
	values[nMonomers*(1+2*nPolymers)+nPolymers:nMonomers*(1+2*nPolymers)+2*nPolymers] = -1
	values[nMonomers*(1+2*nPolymers)+2*nPolymers:nMonomers*(1+2*nPolymers)+4*nPolymers] = 1
	values[nMonomers*(1+2*nPolymers)+4*nPolymers:nMonomers*(1+2*nPolymers)+5*nPolymers] = -1

	rowIndexes[:nMonomers*(1+nPolymers)] = np.tile(np.arange(nMonomers), 1+nPolymers)
	rowIndexes[nMonomers*(1+nPolymers):nMonomers*(1+2*nPolymers)] = nMonomers + np.repeat(np.arange(nPolymers), nMonomers)
	rowIndexes[nMonomers*(1+2*nPolymers):nMonomers*(1+2*nPolymers)+2*nPolymers] = np.tile(nMonomers + np.arange(nPolymers), 2)
	rowIndexes[nMonomers*(1+2*nPolymers)+2*nPolymers:nMonomers*(1+2*nPolymers)+5*nPolymers] = np.tile(nMonomers + nPolymers + np.arange(nPolymers), 3)

	colIndexes[:nMonomers*(1+nPolymers)] = np.arange(nMonomers*(1+nPolymers))
	colIndexes[nMonomers*(1+nPolymers):nMonomers*(1+2*nPolymers)] = nMonomers + np.arange(nMonomers*nPolymers)
	colIndexes[nMonomers*(1+2*nPolymers):nMonomers*(1+2*nPolymers)+nPolymers] = nMonomers*(1+nPolymers) + np.arange(nPolymers)
	colIndexes[nMonomers*(1+2*nPolymers)+nPolymers:nMonomers*(1+2*nPolymers)+2*nPolymers] = nMonomers*(1+nPolymers) + nPolymers
	colIndexes[nMonomers*(1+2*nPolymers)+2*nPolymers:nMonomers*(1+2*nPolymers)+3*nPolymers] = nMonomers*(1+nPolymers) + np.arange(nPolymers)
	colIndexes[nMonomers*(1+2*nPolymers)+3*nPolymers:nMonomers*(1+2*nPolymers)+4*nPolymers] = nMonomers*(1+nPolymers) + nPolymers
	colIndexes[nMonomers*(1+2*nPolymers)+4*nPolymers:nMonomers*(1+2*nPolymers)+5*nPolymers] = nMonomers*(1+nPolymers) + nPolymers + 1 + np.arange(nPolymers)

	temp_matrix = np.zeros((nNodes, nEdges))

	temp_matrix[rowIndexes, colIndexes] = values

	# Build the bounds

	lowerBounds = np.zeros(nEdges)

	upperBounds = np.empty(nEdges)

	upperBounds[:nMonomers] = monomerCounts
	upperBounds[nMonomers:nMonomers+nMonomers*nPolymers] = monomerDeficits.reshape(-1)
	upperBounds[nMonomers+nMonomers*nPolymers:nMonomers*(1+nPolymers)+nPolymers+1] = np.inf
	upperBounds[-nPolymers:] = maxElongation

	# Build the objective

	# Add some noise to prevent systematic bias in polymer NT assignment
	objectiveNoise = 0.01 * randStream.rand(nPolymers)

	# Provide a bonus to polymers that elongate past the global elongation minimum
	bonusElongationWeight = 0.1/nPolymers

	objective = np.zeros(nEdges)
	objective[nMonomers*(1+nPolymers):nMonomers*(1+nPolymers)+nPolymers] = -bonusElongationWeight * (1 + objectiveNoise)
	objective[nMonomers*(1+nPolymers)+nPolymers] = -1

	# Build cvxopt-typed matrices/vectors

	A = spmatrix(values, rowIndexes, colIndexes)
	b = matrix(np.zeros(nNodes))

	G = spmatrix(
		[1]*nEdges + [-1]*nEdges,
		np.arange(2*nEdges),
		np.tile(np.arange(nEdges), 2)
		)

	h = matrix(np.concatenate([upperBounds, -lowerBounds], axis = 0))

	f = matrix(objective)

	# Run (I)LP solver

	if useIntegerLinearProgramming:
		I = set(xrange(nEdges))
		solution = ilp(f, G, h, A, b, I)

	else:
		solution = lp(f, G, h, A, b)

	# Extract solution (monomers assigned)

	v = np.array(solution[1]).astype(np.int).flatten()

	assignments = v[nMonomers:nMonomers+nMonomers*nPolymers]

	monomerAssignments = assignments.reshape(nPolymers, nMonomers)

	return monomerAssignments

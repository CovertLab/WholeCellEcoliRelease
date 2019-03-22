"""
_complexation.pyx

Forms subunits into complexes using random complexation reaction sampling.

TODO:
- document algorithm (not terribly complicated)

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/17/14
"""

from __future__ import absolute_import
from __future__ import division

import numpy as np
cimport numpy as np

from numpy.random import RandomState

np.import_array()

cimport cython

cdef np.int64_t MAX_ITERATIONS = 10**9
cdef np.int64_t NO_MORE_ENTRIES = -1
DEF RANDOM_BATCH_SIZE = 2048

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef mccBuildMatrices(np.ndarray[np.int64_t, ndim=2] stoichiometricMatrix):

	# Collect matrix size information
	cdef int nMolecules = stoichiometricMatrix.shape[0]
	cdef int nReactions = stoichiometricMatrix.shape[1]

	# Collect subunit information
	cdef np.ndarray[np.int64_t, ndim=2] usesMolecule = (stoichiometricMatrix != 0).astype(np.int64)
	cdef int maxMoleculeTypes = np.max(np.sum(usesMolecule, 0))

	cdef np.ndarray[np.int64_t, ndim=2] moleculeIndexes = np.empty((nReactions, maxMoleculeTypes), np.int64)
	moleculeIndexes.fill(NO_MORE_ENTRIES)

	cdef int reactionIndex, moleculeIndex, nSubunits

	for reactionIndex in range(nReactions):
		nSubunits = 0
		for moleculeIndex in range(nMolecules):
			if usesMolecule[moleculeIndex, reactionIndex]:
				moleculeIndexes[reactionIndex, nSubunits] = moleculeIndex
				nSubunits = nSubunits + 1

	# Find which reactions overlap in molecule usage/production
	cdef int maxReactionOverlap = (np.dot(usesMolecule.T, usesMolecule) != 0).sum(axis = 1).max()

	cdef np.ndarray[np.int64_t, ndim=2] overlappingReactions = np.empty((nReactions, maxReactionOverlap), np.int64)
	overlappingReactions.fill(NO_MORE_ENTRIES)

	cdef int subunitIndex, reactionIndex2, nOverlaps
	cdef np.int64_t doesOverlap

	for reactionIndex in range(nReactions):
		nOverlaps = 0
		for reactionIndex2 in range(nReactions):
			doesOverlap = 0
			for subunitIndex in range(maxMoleculeTypes):
				moleculeIndex = moleculeIndexes[reactionIndex, subunitIndex]

				if moleculeIndex == NO_MORE_ENTRIES:
					break

				if usesMolecule[moleculeIndex, reactionIndex2]:
					doesOverlap = 1
					break

			if doesOverlap:
				overlappingReactions[reactionIndex, nOverlaps] = reactionIndex2
				nOverlaps = nOverlaps + 1

	return (moleculeIndexes, overlappingReactions)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef np.ndarray[np.int64_t, ndim=1] mccFormComplexesWithPrebuiltMatrices(
		np.ndarray[np.int64_t, ndim=1] moleculeCounts,
		unsigned int seed,
		np.ndarray[np.int64_t, ndim=2] stoichiometricMatrix,
		np.ndarray[np.int64_t, ndim=2] moleculeIndexes,
		np.ndarray[np.int64_t, ndim=2] overlappingReactions,
		):

	# Seed a random number instance. Get a batch of samples at a time for speed.
	cdef random_state = RandomState(seed)
	cdef np.float64_t[:] random_batch = None
	cdef int random_batch_position = RANDOM_BATCH_SIZE

	# Copy the molecule counts into a new vector for return
	cdef np.ndarray[np.int64_t, ndim=1] updatedMoleculeCounts = moleculeCounts.copy()

	# Collect matrix size information
	cdef int nReactions = stoichiometricMatrix.shape[1]
	cdef int maxMoleculeTypes = moleculeIndexes.shape[1]
	cdef int maxReactionOverlap = overlappingReactions.shape[1]

	# Create vectors for choosing reactions
	cdef np.ndarray[np.int64_t, ndim=1] reactionIsPossible = np.empty(nReactions, np.int64)
	cdef np.ndarray[np.int64_t, ndim=1] reactionCumulative = np.empty_like(reactionIsPossible)

	# Find which reactions are possible
	cdef int reactionIndex, subunitIndex, moleculeIndex
	cdef np.int64_t reactionPossible

	for reactionIndex in range(nReactions):
		reactionPossible = 1
		for subunitIndex in range(maxMoleculeTypes):
			moleculeIndex = moleculeIndexes[reactionIndex, subunitIndex]

			if moleculeIndex == NO_MORE_ENTRIES:
				break

			if updatedMoleculeCounts[moleculeIndex] < -stoichiometricMatrix[moleculeIndex, reactionIndex]:
				reactionPossible = 0
				break

		reactionIsPossible[reactionIndex] = reactionPossible

	# Recursively form complexes
	cdef np.float64_t random_double, cutoffValue, maximumValue
	cdef int iteration, overlapIndex, reactionIndex2

	for iteration in range(MAX_ITERATIONS):

		# Find the total reaction potential
		for reactionIndex in range(nReactions):
			if reactionIndex == 0:
				reactionCumulative[reactionIndex] = reactionIsPossible[reactionIndex]

			else:
				reactionCumulative[reactionIndex] = (
					reactionCumulative[reactionIndex-1]
					+ reactionIsPossible[reactionIndex]
					)

		maximumValue = reactionCumulative[nReactions-1]

		# Break if no reactions are possible
		if maximumValue == 0:
			break

		# Choose which reaction to perform
		if random_batch_position >= RANDOM_BATCH_SIZE:
			random_batch = random_state.random_sample(RANDOM_BATCH_SIZE)
			random_batch_position = 0
		random_double = random_batch[random_batch_position]
		random_batch_position += 1
		cutoffValue = random_double * maximumValue

		for reactionIndex in range(nReactions):
			if cutoffValue <= reactionCumulative[reactionIndex]:
				break

		# Perform the reaction
		for subunitIndex in range(maxMoleculeTypes):
			moleculeIndex = moleculeIndexes[reactionIndex, subunitIndex]

			if moleculeIndex == NO_MORE_ENTRIES:
				break

			updatedMoleculeCounts[moleculeIndex] = (
				updatedMoleculeCounts[moleculeIndex] +
				stoichiometricMatrix[moleculeIndex, reactionIndex]
				)

		# Update relevant reactions
		for overlapIndex in range(maxReactionOverlap):
			reactionIndex2 = overlappingReactions[reactionIndex, overlapIndex]

			if reactionIndex2 == NO_MORE_ENTRIES:
				break

			reactionPossible = 1
			for subunitIndex in range(maxMoleculeTypes):
				moleculeIndex = moleculeIndexes[reactionIndex2, subunitIndex]

				if moleculeIndex == NO_MORE_ENTRIES:
					break

				if updatedMoleculeCounts[moleculeIndex] < -stoichiometricMatrix[moleculeIndex, reactionIndex2]:
					reactionPossible = 0
					break

			reactionIsPossible[reactionIndex2] = reactionPossible

	return updatedMoleculeCounts


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef np.ndarray[np.int64_t, ndim=1] mccFormComplexes(
		np.ndarray[np.int64_t, ndim=1] moleculeCounts,
		int seed,
		np.ndarray[np.int64_t, ndim=2] stoichiometricMatrix,
		):

	cdef np.ndarray[np.int64_t, ndim=2] moleculeIndexes, overlappingReactions
	cdef np.ndarray[np.int64_t, ndim=1] updatedMoleculeCounts

	moleculeIndexes, overlappingReactions = mccBuildMatrices(stoichiometricMatrix)

	updatedMoleculeCounts = mccFormComplexesWithPrebuiltMatrices(
		moleculeCounts,
		seed,
		stoichiometricMatrix,
		moleculeIndexes,
		overlappingReactions
		)

	return updatedMoleculeCounts

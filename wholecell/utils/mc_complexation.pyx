"""
_complexation.pyx

Forms subunits into complexes using random complexation reaction sampling.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/17/14
"""

from __future__ import division

import numpy as np
cimport numpy as np

from libc.stdlib cimport rand, srand, RAND_MAX

np.import_array()

cimport cython

cdef np.int64_t MAX_ITERATIONS = 10**9
cdef np.int64_t NO_MORE_ENTRIES = -1

# needed for iteration:
# stoichiometricMatrix, nReactions
# maxSubunitTypes, subunitIndexes
# maxReactionOverlap, overlappingReactions

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef mccBuildMatrices(np.ndarray[np.int64_t, ndim=2] stoichiometricMatrix):

	# Collect matrix size information
	cdef int nMolecules = stoichiometricMatrix.shape[0]
	cdef int nReactions = stoichiometricMatrix.shape[1]

	# Collect subunit information
	cdef np.ndarray[np.int64_t, ndim=2] usesSubunit = (stoichiometricMatrix < 0).astype(np.int64)
	cdef int maxSubunitTypes = np.max(np.sum(usesSubunit, 0))

	cdef np.ndarray[np.int64_t, ndim=2] subunitIndexes = np.empty((nReactions, maxSubunitTypes), np.int64)
	subunitIndexes.fill(NO_MORE_ENTRIES)

	cdef int reactionIndex, moleculeIndex, nSubunits

	for reactionIndex in range(nReactions):
		nSubunits = 0
		for moleculeIndex in range(nMolecules):
			if usesSubunit[moleculeIndex, reactionIndex]:
				subunitIndexes[reactionIndex, nSubunits] = moleculeIndex
				nSubunits = nSubunits + 1

	# Find which reactions overlap in molecule usage/production

	cdef np.ndarray[np.int64_t, ndim=2] usesMolecule = (stoichiometricMatrix != 0).astype(np.int64)
	cdef int maxReactionOverlap = np.max(np.sum(usesMolecule, 0))

	cdef np.ndarray[np.int64_t, ndim=2] overlappingReactions = np.empty((nReactions, maxReactionOverlap), np.int64)
	overlappingReactions.fill(NO_MORE_ENTRIES)

	cdef int subunitIndex, reactionIndex2, nOverlaps
	cdef np.int64_t doesOverlap

	for reactionIndex in range(nReactions):
		nOverlaps = 0
		for reactionIndex2 in range(nReactions):
			doesOverlap = 0
			for subunitIndex in range(maxSubunitTypes):
				moleculeIndex = subunitIndexes[reactionIndex, subunitIndex]

				if moleculeIndex == NO_MORE_ENTRIES:
					break

				if usesMolecule[moleculeIndex, reactionIndex2]:
					doesOverlap = 1
					break

			if doesOverlap:
				overlappingReactions[reactionIndex, nOverlaps] = reactionIndex2
				nOverlaps = nOverlaps + 1

	return (subunitIndexes, overlappingReactions)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef np.ndarray[np.int64_t, ndim=1] mccFormComplexesWithPrebuiltMatrices(
		np.ndarray[np.int64_t, ndim=1] moleculeCounts,
		int seed,
		np.ndarray[np.int64_t, ndim=2] stoichiometricMatrix,
		np.ndarray[np.int64_t, ndim=2] subunitIndexes,
		np.ndarray[np.int64_t, ndim=2] overlappingReactions,
		):

	# Set the seed
	srand(seed)
	
	# Copy the molecule counts into a new vector for return
	cdef np.ndarray[np.int64_t, ndim=1] updatedMoleculeCounts = moleculeCounts.copy()

	# Collect matrix size information
	cdef int nReactions = stoichiometricMatrix.shape[1]
	cdef int maxSubunitTypes = subunitIndexes.shape[1]
	cdef int maxReactionOverlap = overlappingReactions.shape[1]

	# Create vectors for choosing reactions
	cdef np.ndarray[np.int64_t, ndim=1] reactionIsPossible = np.empty(nReactions, np.int64)
	cdef np.ndarray[np.int64_t, ndim=1] reactionCumulative = np.empty_like(reactionIsPossible)

	# Find which reactions are possible
	cdef int reactionIndex, subunitIndex, moleculeIndex
	cdef np.int64_t reactionPossible

	for reactionIndex in range(nReactions):
		reactionPossible = 1
		for subunitIndex in range(maxSubunitTypes):
			moleculeIndex = subunitIndexes[reactionIndex, subunitIndex]

			if moleculeIndex == NO_MORE_ENTRIES:
				break

			if updatedMoleculeCounts[moleculeIndex] < -stoichiometricMatrix[moleculeIndex, reactionIndex]:
				reactionPossible = 0
				break

		reactionIsPossible[reactionIndex] = reactionPossible

	# Recursively form complexes
	cdef np.float64_t random, cutoffValue, maximumValue
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
		random = <np.float64_t>rand() / <np.float64_t>RAND_MAX
		cutoffValue = random * maximumValue

		for reactionIndex in range(nReactions):
			if cutoffValue <= reactionCumulative[reactionIndex]:
				break

		# Perform the reaction
		for subunitIndex in range(maxSubunitTypes):
			moleculeIndex = subunitIndexes[reactionIndex, subunitIndex]

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
			for subunitIndex in range(maxSubunitTypes):
				moleculeIndex = subunitIndexes[reactionIndex2, subunitIndex]

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
		np.ndarray[np.int64_t, ndim=2] stoichiometricMatrix,
		int seed,
		):

	cdef np.ndarray[np.int64_t, ndim=2] subunitIndexes, overlappingReactions
	cdef np.ndarray[np.int64_t, ndim=1] updatedMoleculeCounts

	subunitIndexes, overlappingReactions = mccBuildMatrices(stoichiometricMatrix)

	updatedMoleculeCounts = mccFormComplexesWithPrebuiltMatrices(
		moleculeCounts,
		seed,
		stoichiometricMatrix,
		subunitIndexes,
		overlappingReactions
		)

	return updatedMoleculeCounts

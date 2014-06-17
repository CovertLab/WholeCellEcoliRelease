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

cdef np.int64_t DEFAULT_MAX_ITERATIONS = 10**9
cdef np.int64_t NO_MORE_ENTRIES = -1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef np.ndarray[np.int64_t, ndim=1] monteCarloComplexation(
		np.ndarray[np.int64_t, ndim=1] moleculeCounts,
		np.ndarray[np.int64_t, ndim=2] stoichiometricMatrix,
		int seed,
		np.int64_t maxIterations = DEFAULT_MAX_ITERATIONS,
		):

	# Set the seed
	srand(seed)
	
	# Copy the molecule counts into a new vector for return
	cdef np.ndarray[np.int64_t, ndim=1] updatedMoleculeCounts = moleculeCounts.copy()

	# Collect matrix size information
	cdef np.int64_t nMolecules = stoichiometricMatrix.shape[0]
	cdef np.int64_t nReactions = stoichiometricMatrix.shape[1]

	# Collect subunit information
	cdef np.ndarray[np.int64_t, ndim=2] usesSubunit = (stoichiometricMatrix < 0).astype(np.int64)
	cdef np.int64_t maxSubunitTypes = np.max(np.sum(usesSubunit, 0))

	cdef np.ndarray[np.int64_t, ndim=2] subunitIndexes = np.empty((nReactions, maxSubunitTypes), np.int64)
	subunitIndexes.fill(NO_MORE_ENTRIES)

	cdef int reactionIndex, moleculeIndex, nSubunits

	for reactionIndex in range(nReactions):
		nSubunits = 0
		for moleculeIndex in range(nMolecules):
			if usesSubunit[moleculeIndex, reactionIndex]:
				subunitIndexes[reactionIndex, nSubunits] = moleculeIndex
				nSubunits = nSubunits + 1

	del usesSubunit

	# Find which reactions overlap in molecule usage/production

	cdef np.ndarray[np.int64_t, ndim=2] usesMolecule = (stoichiometricMatrix != 0).astype(np.int64)
	cdef np.int64_t maxReactionOverlap = np.max(np.sum(usesMolecule, 0))

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

				assert nOverlaps < maxReactionOverlap

	del usesMolecule

	# Create vectors for choosing reactions
	cdef np.ndarray[np.int64_t, ndim=1] reactionIsPossible = np.empty(nReactions, np.int64)
	cdef np.ndarray[np.int64_t, ndim=1] reactionCumulative = np.empty_like(reactionIsPossible)

	# Find which reactions are possible
	cdef int reactionPossible

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
	cdef int iteration, overlapIndex

	for iteration in range(maxIterations):

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

		# else:
		# 	raise Exception("Failed to find a valid reaction - this should never happen")

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

		# if <int>(<float>iteration % 10000.) == 0:
		# 	print iteration, np.sum(updatedMoleculeCounts), reactionIndex, maximumValue
		# 	assert (updatedMoleculeCounts >= 0).all()

	# else:
	# 	raise Exception("Reached maxium iteration limit")

	return updatedMoleculeCounts

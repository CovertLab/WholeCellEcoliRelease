# cython: language_level=3str
## [Enable this in Cython 3]  distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

"""
_complexation.pyx

Forms subunits into complexes using random complexation reaction sampling.

TODO:
- document algorithm (not terribly complicated)
"""

from __future__ import absolute_import, division, print_function

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
	cdef const np.int64_t[:, :] usesMoleculeView = usesMolecule  # it's not always C-contiguous
	cdef int maxMoleculeTypes = np.max(np.sum(usesMolecule, 0))

	cdef np.ndarray[np.int64_t, ndim=2] moleculeIndexes = np.full(
		(nReactions, maxMoleculeTypes), NO_MORE_ENTRIES, np.int64)

	cdef int reactionIndex, moleculeIndex, nSubunits

	for reactionIndex in range(nReactions):
		nSubunits = 0
		for moleculeIndex in range(nMolecules):
			if usesMoleculeView[moleculeIndex, reactionIndex]:
				moleculeIndexes[reactionIndex, nSubunits] = moleculeIndex
				nSubunits = nSubunits + 1

	# Find which reactions overlap in molecule usage/production
	cdef int maxReactionOverlap = (np.dot(usesMolecule.T, usesMolecule) != 0).sum(axis = 1).max()

	cdef np.ndarray[np.int64_t, ndim=2] overlappingReactions = np.full(
		(nReactions, maxReactionOverlap), NO_MORE_ENTRIES, np.int64)
	cdef np.int64_t[:, ::1] overlappingReactionsView = overlappingReactions

	cdef int subunitIndex, reactionIndex2, nOverlaps, doesOverlap

	for reactionIndex in range(nReactions):
		nOverlaps = 0
		for reactionIndex2 in range(nReactions):
			doesOverlap = 0
			for subunitIndex in range(maxMoleculeTypes):
				moleculeIndex = moleculeIndexes[reactionIndex, subunitIndex]

				if moleculeIndex == NO_MORE_ENTRIES:
					break

				if usesMoleculeView[moleculeIndex, reactionIndex2]:
					doesOverlap = 1
					break

			if doesOverlap:
				overlappingReactionsView[reactionIndex, nOverlaps] = reactionIndex2
				nOverlaps += 1

	return moleculeIndexes, overlappingReactions


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef tuple mccFormComplexesWithPrebuiltMatrices(
		np.ndarray[np.int64_t, ndim=1] moleculeCounts,
		unsigned int seed,
		np.ndarray[np.int64_t, ndim=2] stoichiometricMatrix,
		np.ndarray[np.int64_t, ndim=2] moleculeIndexes,
		np.ndarray[np.int64_t, ndim=2] overlappingReactions):

	# Seed a random number instance. Get a batch of samples at a time for speed.
	cdef random_state = RandomState(seed)
	cdef np.float64_t[:] random_batch = None
	cdef int random_batch_position = RANDOM_BATCH_SIZE

	# Copy the molecule counts into a new vector for return
	cdef np.ndarray[np.int64_t, ndim=1] updatedMoleculeCounts = moleculeCounts.copy()
	cdef np.int64_t[::1] updatedMoleculeCountsView = updatedMoleculeCounts

	# Collect matrix size information
	cdef int nReactions = stoichiometricMatrix.shape[1]
	cdef int maxMoleculeTypes = moleculeIndexes.shape[1]
	cdef int maxReactionOverlap = overlappingReactions.shape[1]

	# Create vectors for choosing reactions
	cdef np.ndarray[np.int64_t, ndim=1] reactionIsPossible = np.empty(nReactions, np.int64)
	cdef np.int64_t[::1] reactionIsPossibleView = reactionIsPossible
	cdef np.int64_t[::1] reactionCumulative = np.empty_like(reactionIsPossible)

	# Create vector for saving number of events for each complexation reaction (number of reactions / time step)
	cdef np.ndarray[np.int64_t, ndim=1] complexationEvents = np.zeros(nReactions, np.int64)

	# Find which reactions are possible
	cdef int reactionIndex, subunitIndex, moleculeIndex
	cdef np.int64_t reactionPossible

	for reactionIndex in range(nReactions):
		reactionPossible = 1
		for subunitIndex in range(maxMoleculeTypes):
			moleculeIndex = moleculeIndexes[reactionIndex, subunitIndex]

			if moleculeIndex == NO_MORE_ENTRIES:
				break

			if updatedMoleculeCountsView[moleculeIndex] < -stoichiometricMatrix[moleculeIndex, reactionIndex]:
				reactionPossible = 0
				break

		reactionIsPossibleView[reactionIndex] = reactionPossible

	# Recursively form complexes
	cdef np.float64_t random_double, cutoffValue, maximumValue
	cdef int iteration, overlapIndex, reactionIndex2

	for iteration in range(MAX_ITERATIONS):

		# Find the total reaction potential
		for reactionIndex in range(nReactions):
			if reactionIndex == 0:
				reactionCumulative[reactionIndex] = reactionIsPossibleView[reactionIndex]

			else:
				reactionCumulative[reactionIndex] = (
					reactionCumulative[reactionIndex-1] + reactionIsPossibleView[reactionIndex])

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

			updatedMoleculeCountsView[moleculeIndex] = (
				updatedMoleculeCountsView[moleculeIndex] +
				stoichiometricMatrix[moleculeIndex, reactionIndex])

			complexationEvents[reactionIndex] += 1

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

				if updatedMoleculeCountsView[moleculeIndex] < -stoichiometricMatrix[moleculeIndex, reactionIndex2]:
					reactionPossible = 0
					break

			reactionIsPossibleView[reactionIndex2] = reactionPossible

	return updatedMoleculeCounts, complexationEvents


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef np.ndarray[np.int64_t, ndim=1] mccFormComplexes(
		np.ndarray[np.int64_t, ndim=1] moleculeCounts,
		unsigned int seed,
		np.ndarray[np.int64_t, ndim=2] stoichiometricMatrix):

	cdef np.ndarray[np.int64_t, ndim=2] moleculeIndexes, overlappingReactions
	cdef np.ndarray[np.int64_t, ndim=1] updatedMoleculeCounts, complexationEvents

	moleculeIndexes, overlappingReactions = mccBuildMatrices(stoichiometricMatrix)

	updatedMoleculeCounts, complexationEvents = mccFormComplexesWithPrebuiltMatrices(
		moleculeCounts,
		seed,
		stoichiometricMatrix,
		moleculeIndexes,
		overlappingReactions
		)

	return updatedMoleculeCounts

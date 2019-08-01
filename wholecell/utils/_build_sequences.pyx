"""
_build_sequences.pyx

Builds the matrices used for the polymerize function.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/23/14
"""

from __future__ import division

import numpy as np

import numpy as np
cimport numpy as np

np.import_array()

cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef np.ndarray[np.int8_t, ndim=2] buildSequences(
		np.ndarray[np.int8_t, ndim=2] base_sequences,
		np.ndarray[np.int64_t, ndim=1] indexes,
		np.ndarray[np.int64_t, ndim=1] positions,
		np.ndarray[np.int64_t, ndim=1] elongation_rates):

	cdef int PAD_VALUE = -1

	cdef int elongation_max = elongation_rates.max()
	if np.any(positions + elongation_max > base_sequences.shape[1]):
		raise Exception('Elongation proceeds past end of sequence!')

	cdef int out_rows = positions.shape[0]

	cdef np.ndarray[np.int8_t, ndim=2] out = np.empty((out_rows, elongation_max), np.int8)

	cdef int i, index, position, j

	for i in range(out_rows):
		index = indexes[i]
		position = positions[i]
		rate = elongation_rates[index]

		for j in range(elongation_max):
			if j < rate:
				out[i, j] = base_sequences[index, position+j]
			else:
				out[i, j] = PAD_VALUE

	return out

# TODO: move this to a new file?

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef np.ndarray[np.float64_t, ndim=1] computeMassIncrease(
		np.ndarray[np.int8_t, ndim=2] sequences,
		np.ndarray[np.int64_t, ndim=1] elongations,
		np.ndarray[np.float64_t, ndim=1] monomerMasses):

	cdef int out_size = sequences.shape[0]

	cdef np.ndarray[np.float64_t, ndim=1] out = np.zeros(out_size, np.float64)

	cdef int i, j

	for i in range(out_size):
		for j in range(elongations[i]):
			out[i] += monomerMasses[
				sequences[i, j]]

	return out

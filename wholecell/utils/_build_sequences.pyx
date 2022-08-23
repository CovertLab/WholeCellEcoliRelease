# cython: language_level=3str
## [Enable this in Cython 3]  distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

"""
_build_sequences.pyx

Builds the matrices used for the polymerize function.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
cimport numpy as np

np.import_array()

cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef np.ndarray[np.int8_t, ndim=2] buildSequences(
		const np.int8_t[:, ::1] base_sequences,
		const np.int64_t[::1] indexes,
		np.ndarray[np.int64_t, ndim=1] positions,
		np.ndarray[np.int64_t, ndim=1] elongation_rates):

	cdef int elongation_max = elongation_rates.max()
	if np.any(positions + elongation_max > base_sequences.shape[1]):
		raise IndexError('Elongation proceeds past end of sequence!')

	cdef const np.int64_t[::1] positions_view = positions
	cdef int out_rows = positions_view.shape[0]

	cdef np.ndarray[np.int8_t, ndim=2] out = np.empty((out_rows, elongation_max), np.int8)
	cdef np.int8_t[:, ::1] out_view = out

	cdef int i, index, position, j

	for i in range(out_rows):
		index = indexes[i]
		position = positions_view[i]

		for j in range(elongation_max):
			out_view[i, j] = base_sequences[index, position+j]

	return out

# TODO: move this to a new file?

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef np.ndarray[np.float64_t, ndim=1] computeMassIncrease(
		const np.int8_t[:, ::1] sequences,
		const np.int64_t[::1] elongations,
		const np.float64_t[::1] monomerMasses):

	cdef int out_size = sequences.shape[0]

	cdef np.ndarray[np.float64_t, ndim=1] out = np.zeros(out_size, np.float64)
	cdef np.float64_t[::1] out_view = out

	cdef int i, j

	for i in range(out_size):
		for j in range(elongations[i]):
			out_view[i] = out_view[i] + monomerMasses[sequences[i, j]]

	return out

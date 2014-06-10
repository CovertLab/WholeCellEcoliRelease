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

DTYPE8 = np.int8
ctypedef np.int8_t DTYPE8_t
ctypedef np.int64_t DTYPE64_t

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef np.ndarray[DTYPE8_t, ndim=2] buildSequences(
		np.ndarray[DTYPE8_t, ndim=2] base_sequences,
		np.ndarray[DTYPE64_t, ndim=1] indexes,
		np.ndarray[DTYPE64_t, ndim=1] positions,
		int elng_rate):

	cdef int out_rows = positions.shape[0]

	cdef np.ndarray[DTYPE8_t, ndim=2] out = np.empty((out_rows, elng_rate), DTYPE8)

	cdef int i, index, position, j

	for i in range(out_rows):
		index = indexes[i]
		position = positions[i]

		for j in range(elng_rate):
			out[i, j] = base_sequences[index, position+j]

	return out

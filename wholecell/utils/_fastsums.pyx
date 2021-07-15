# cython: language_level=3str
## [Enable this in Cython 3]  distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

"""
_fastsums.pyx

Fast group sums to speed up polymerize.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
cimport numpy as np
cimport cython


ctypedef Py_ssize_t Index # array index type
ctypedef np.int32_t Int32

def sum_monomers_reference_implementation(sequenceMonomers, activeSequencesIndexes):
	"""
	Sum up the total number of monomers of each type needed to continue building
	the active sequences through currentStep. (This is the
	Python reference implementation, compiled by Cython.)

	Arguments:
	sequenceMonomers -- bool[monomer #, sequence #] indicating whether
		a given monomer gets used in a step of a sequence
	activeSequencesIndexes -- an array of sequences that are still active, i.e.
		have not yet run out of source monomers

	Result:
	count[monomer #] indicating how many of each monomer will be needed
		by the combined active sequences in the currentStep
	"""
	totalMonomers = (
		sequenceMonomers[:, activeSequencesIndexes] # filter to active sequences
		.sum(axis=1)) # sum over all active sequences
	return totalMonomers

def sum_monomers(
		sequenceMonomers not None,
		const np.int64_t[::1] monomerIndexes not None,
		const np.int64_t[::1] activeSequencesIndexes not None):
	"""
	Sum up the total number of monomers of each type needed to continue building
	the active sequences for the currentStep.

	Workaround: Cython doesn't support boolean array Memoryviews, so accept an
	untyped sequenceMonomers array and pass a uint8[] view of it to _sum_monomers().
	"""
	return _sum_monomers(
		sequenceMonomers.view(dtype=np.uint8),
		monomerIndexes,
		activeSequencesIndexes)

# See:
# http://cython.readthedocs.io/en/latest/src/tutorial/numpy.html
# http://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html
# http://www.perrygeo.com/parallelizing-numpy-array-loops-with-cython-and-mpi.html

# TODO:
#   Make activeSequencesIndexes an int32[] instead of int64[].
#   http://cython.readthedocs.io/en/latest/src/reference/compilation.html#compiler-directives
#   Try transposing sumMonomer axis 1 & 2 for better memory/cache locality.
#   Pack 8 booleans into a byte for better memory/cache performance.
#   Release the gil and parallelize the outer loop.
@cython.wraparound(False)   # --> Without boundscheck(False) this is slower!
@cython.boundscheck(False)
cdef np.ndarray _sum_monomers(
		const np.uint8_t[:, :, ::1] sequenceMonomers,
		const np.int64_t[::1] monomerIndexes,
		const np.int64_t[::1] activeSequencesIndexes):
	cdef Index nMonomers = sequenceMonomers.shape[0]
	cdef Index maxSequences = sequenceMonomers.shape[1]
	cdef Index nActiveSequences = activeSequencesIndexes.shape[0]
	cdef Int32 total = 0
	cdef Index monomer, step, iseq, seq

	# Do the bounds-checks once before looping since @cython.boundscheck(False).
	for iseq in range(nActiveSequences):
		seq = activeSequencesIndexes[iseq]
		if seq < 0 or seq >= maxSequences:
			raise IndexError('activeSequencesIndexes[%s]=%s is out of range(%s)'
				% (iseq, seq, maxSequences))

	cdef np.ndarray totalMonomers = np.empty(nMonomers, dtype=np.int32)
	cdef Int32[:] _totalMonomers = totalMonomers # a typed memoryview of totalMonomers

	for monomer in range(nMonomers):
		total = 0
		for iseq in range(nActiveSequences):
			seq = activeSequencesIndexes[iseq]
			total += sequenceMonomers[monomer, seq, monomerIndexes[iseq]]
		_totalMonomers[monomer] = total

	return totalMonomers

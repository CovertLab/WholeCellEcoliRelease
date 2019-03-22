"""
_fastsums.pyx

Fast group sums to speed up polymerize.

@author: Jerry Morrison
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 10/12/2016
"""

from __future__ import division

import numpy as np
cimport numpy as np
cimport cython


ctypedef Py_ssize_t Index # array index type
ctypedef np.int32_t Int32

def sum_monomers_reference_implementation(
		sequenceMonomers, activeSequencesIndexes, Index currentStep):
	"""
	Sum up the total number of monomers of each type needed to continue building
	the active sequences through currentStep and following steps. (This is the
	Python reference implementation, compiled by Cython.)

	Arguments:
	sequenceMonomers -- bool[monomer #, sequence #, step #] indicating whether
		a given monomer gets used in a step of a sequence
	activeSequencesIndexes -- an array of sequences that are still active, i.e.
		have not yet run out of source monomers
	currentStep -- the current sequence step number

	Result:
	count[monomer #, step #] indicating how many of each monomer will be needed
		by the combined active sequences in the currentStep and following steps
	"""
	totalMonomers = (
		sequenceMonomers[:, activeSequencesIndexes, currentStep:] # filter to active sequences
		.sum(axis=1) # sum over all active sequences
		.cumsum(axis=1)) # cumsum from currentStep to the last step
	return totalMonomers

def sum_monomers(
		sequenceMonomers, np.int64_t[:] activeSequencesIndexes, currentStep):
	"""
	Sum up the total number of monomers of each type needed to continue building
	the active sequences through currentStep and following steps.

	Workaround: Cython doesn't support boolean array Memoryviews, so accept an
	untyped sequenceMonomers array and pass a uint8[] view of it to _sum_monomers().
	"""
	return _sum_monomers(
		sequenceMonomers.view(dtype=np.uint8), activeSequencesIndexes, currentStep)

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
def _sum_monomers(
		np.uint8_t[:, :, ::1] sequenceMonomers not None,
		np.int64_t[::1] activeSequencesIndexes not None,
		Index currentStep):
	cdef Index nMonomers = sequenceMonomers.shape[0]
	cdef Index maxSequences = sequenceMonomers.shape[1]
	cdef Index maxSteps = sequenceMonomers.shape[2]
	cdef Index nActiveSequences = activeSequencesIndexes.shape[0]
	cdef Int32 total = 0
	cdef Index monomer, step, iseq, seq

	# Do the bounds-checks once before looping *if* @cython.boundscheck(False).
	# Note: Testing 'currentStep not in xrange(maxSteps)' is very slow!
	if currentStep < 0 or currentStep >= maxSteps:
		raise IndexError('currentStep %s is out of range(%s)' % (currentStep, maxSteps))
	for iseq in xrange(nActiveSequences):
		seq = activeSequencesIndexes[iseq]
		if seq < 0 or seq >= maxSequences:
			raise IndexError('activeSequencesIndexes[%s]=%s is out of range(%s)'
				% (iseq, seq, maxSequences))

	cdef np.ndarray totalMonomers = np.empty((nMonomers, maxSteps - currentStep), dtype=np.int32)
	cdef Int32[:, :] _totalMonomers = totalMonomers # a typed memoryview of totalMonomers

	for monomer in xrange(nMonomers):
		total = 0
		for step in xrange(currentStep, maxSteps):
			for iseq in xrange(nActiveSequences):
				seq = activeSequencesIndexes[iseq]
				total += sequenceMonomers[monomer, seq, step]
			_totalMonomers[monomer, step - currentStep] = total

	return totalMonomers

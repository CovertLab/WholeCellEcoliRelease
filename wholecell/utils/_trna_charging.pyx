# cython: language_level=3str
## [Enable this in Cython 3]  distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

"""
_trna_charging.pyx

These programs are used to speed up steps in KineticTrnaChargingModel.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
cimport numpy as np

from libc.stdlib cimport rand, srand, labs, malloc, free, realloc
from libc.math cimport floor, ceil
from libc.time cimport time, time_t

np.import_array()
cimport cython

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
@cython.cdivision(True)		# Deactivate 0 division checking
cpdef int get_initiations(
	long [:] elongations,
	long [:] lengths,
	long [:] indexes,
	):

	cdef int n_initiations = 0

	# Py_ssize_t is the proper C type for Python array indices.
	cdef Py_ssize_t size = elongations.shape[0]
	cdef Py_ssize_t i

	for i in range(size):
		if elongations[i] > 0 and lengths[i] == 0:
			n_initiations += 1

	return n_initiations

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
@cython.cdivision(True)		# Deactivate 0 division checking
cpdef char get_codon_at(
	char [:, :] sequences, # 8-bit signed integer
	long [:] elongations,
	long ith_ribosome,
	long relative_position,
	int absolute_position,
	):

	absolute_position = elongations[ith_ribosome] - 1 + relative_position

	if absolute_position < 0:
		return -1
	elif absolute_position >= sequences.shape[1]:
		return -1
	else:
		return sequences[ith_ribosome, absolute_position]

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
@cython.cdivision(True)		# Deactivate 0 division checking
cpdef (int, int) get_candidates_to_C(
	char [:, :] sequences,
	long [:] elongations,
	char codon_id,
	int candidates,
	int relative_position,
	int i,
	int absolute_position,
	):

	candidates = 0

	for relative_position in range(1, sequences.shape[1] + 1, 1):
		for i in range(sequences.shape[0]):
			if get_codon_at(
					sequences, elongations, i, relative_position,
					absolute_position) == codon_id:
				candidates += 1

		if candidates > 0:
			break

	return candidates, relative_position


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
@cython.cdivision(True)		# Deactivate 0 division checking
cpdef (int, int) get_candidates_to_N(
	char [:, :] sequences,
	long [:] elongations,
	char codon_id,
	int candidates,
	int relative_position,
	int i,
	int absolute_position,
	):

	candidates = 0

	for relative_position in range(0, -sequences.shape[1], -1):
		for i in range(sequences.shape[0]):
			if get_codon_at(
					sequences, elongations, i, relative_position,
					absolute_position) == codon_id:
				candidates += 1

		if candidates > 0:
			break

	return candidates, relative_position


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
@cython.cdivision(True)		# Deactivate 0 division checking
cpdef int select_candidate(
	char [:, :] sequences,
	long [:] elongations,
	int relative_position,
	char codon_id,
	int r,
	int i,
	int j,
	int absolute_position,
	):

	j = -1
	for i in range(sequences.shape[0]):
		if get_codon_at(
				sequences, elongations, i, relative_position,
				absolute_position) == codon_id:
			j += 1
			if j == r:
				break
	return i


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
@cython.cdivision(True)		# Deactivate 0 division checking
cdef bint is_initial_state(
	int [:] initial_state,
	int [:] state,
	char size,
	int i,
	):
	for i in range(size):
		if initial_state[i] != state[i]:
			return False
	return True


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
@cython.cdivision(True)		# Deactivate 0 division checking
cpdef reconcile_via_ribosome_positions(
	long [:] sequence_codons, # 64-bit integer
	long [:] elongations, # 64-bit integer
	long [:] kinetics_codons, # 64-bit integer
	char [:, :] sequences, # 8-bit integer
	char max_attempts, # 8-bit integer
	):

	# Set random seed
	cdef time_t t = time(NULL)
	srand(t)

	# Initialize variables
	cdef char attempt = 0
	cdef int i = 0
	cdef int j = 0
	cdef int r = 0
	cdef int codon = 0
	cdef int codons = sequence_codons.shape[0]
	cdef int candidates = 0
	cdef int compromise = 0
	cdef int absolute_position = 0
	cdef int relative_position = 0
	cdef int disagreements = 0
	cdef bint disagreements_remaining = True

	# Initialize array of codons that have been exhausted (no longer
	# found in the sequence matrix for this time step)
	cdef char *exhausted = <char *>malloc(codons * sizeof(char))
	if not exhausted:
		raise MemoryError()


	# Search for the best compromise that honors the limits of the
	# Kinetic Model.
	for attempt in range(max_attempts):

		for i in range(codons):
			exhausted[i] = 0

		# Starting from the Sequence Model's solution, take forward
		# steps to match the Kinetic Model's solution.
		while disagreements_remaining:

			# Determine number of disagreements
			disagreements = 0
			for i in range(codons):
				if kinetics_codons[i] > sequence_codons[i] and exhausted[i] == 0:
					disagreements += kinetics_codons[i] - sequence_codons[i]
			if disagreements == 0:
				disagreements_remaining = False
				continue

			# Randomly select a codon disagreement to reconcile
			# Note: prioritizes codon with larger disagreements
			r = rand() % disagreements
			disagreements = -1
			for codon in range(codons):
				if kinetics_codons[codon] > sequence_codons[codon] and exhausted[codon] == 0:
					disagreements += kinetics_codons[codon] - sequence_codons[codon]
					if disagreements >= r:
						break

			# Identify candidate ribosome positions to change
			candidates, relative_position = get_candidates_to_C(
				sequences,
				elongations,
				codon,

				# Variables that are being recycled
				candidates,
				relative_position,
				i,
				absolute_position)
			if candidates == 0:
				exhausted[codon] = 1
				continue

			# Select a candidate
			r = rand() % candidates
			selected_candidate = select_candidate(
				sequences,
				elongations,
				relative_position,
				codon,
				r,

				# Variables that are being recycled
				i,
				j,
				absolute_position)

			# Take forward steps
			for i in range(relative_position):

				codon = get_codon_at(
					sequences,
					elongations,
					selected_candidate,
					1,
					absolute_position)

				# Update Sequence Model
				elongations[selected_candidate] += 1
				sequence_codons[codon] += 1


		# Take backward steps to match the Kinetic Model.
		disagreements_remaining = True
		while disagreements_remaining:

			# Determine number of disagreements
			# Note: No need to monitor exhausted codons because codons
			# must exist in the sequence if the Sequence Model's
			# solution is greater than the Kinetic Model's solution.
			disagreements = 0
			for i in range(codons):
				if kinetics_codons[i] < sequence_codons[i]:
					disagreements += sequence_codons[i] - kinetics_codons[i]
			if disagreements == 0:
				disagreements_remaining = False
				continue

			# Randomly select a codon disagreement to reconcile
			# Note: prioritizes codon with larger disagreements
			r = rand() % disagreements
			disagreements = -1
			for codon in range(codons):
				if kinetics_codons[codon] < sequence_codons[codon]:
					disagreements += sequence_codons[codon] - kinetics_codons[codon]
					if disagreements >= r:
						break

			# Identify candidate ribosome positions to change
			candidates, relative_position = get_candidates_to_N(
				sequences,
				elongations,
				codon,

				# Variables that are being recycled
				candidates,
				relative_position,
				i,
				absolute_position)

			# Select a candidate
			r = rand() % candidates
			selected_candidate = select_candidate(
				sequences,
				elongations,
				relative_position,
				codon,
				r,

				# Variables that are being recycled
				i,
				j,
				absolute_position)

			# Take backward steps
			for i in range(1, relative_position, -1):

				codon = get_codon_at(
					sequences,
					elongations,
					selected_candidate,
					0,
					absolute_position)

				# Update Sequence Model
				elongations[selected_candidate] -= 1
				sequence_codons[codon] -= 1


		# Calculate the compromise
		compromise = 0
		for i in range(codons):
			compromise += labs(kinetics_codons[i] - sequence_codons[i])

		# Accept the compromise
		if compromise == 0:
			break

	free(exhausted)
	return


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
@cython.cdivision(True)		# Deactivate 0 division checking
cpdef reconcile_via_trna_pools(
	long [:] sequence_codons, # 64-bit integer
	long [:] kinetics_codons, # 64-bit integer
	long [:] free_trnas, # 64-bit integer
	long [:] charged_trnas, # 64-bit integer
	long [:] chargings, # 64-bit integer
	long [:] amino_acids_used, # 64-bit integer
	long [:, :]codons_to_trnas_counter, # 64-bit integer
	char [:, :] trnas_to_codons, # 8-bit integer
	char [:] trnas_to_amino_acid_indexes, # 8-bit integer
	):

	# Set random seed
	cdef time_t t = time(NULL)
	srand(t)

	# Initialize variables
	cdef int i = 0
	# cdef int j = 0
	cdef int r = 0
	cdef int codon = 0
	cdef int codons = sequence_codons.shape[0]
	cdef int candidates = 0
	cdef int disagreements = 0
	cdef bint disagreements_remaining = True

	# Reconcile disagreements using tRNAs
	while disagreements_remaining:

		# Determine number of disagreements
		disagreements = 0
		for i in range(codons):
			if kinetics_codons[i] > sequence_codons[i]:
				disagreements += kinetics_codons[i] - sequence_codons[i]
		if disagreements == 0:
			disagreements_remaining = False
			continue

		# Randomly select a codon disagreement to reconcile
		# Note: prioritizes codon with larger disagreements
		r = rand() % disagreements
		disagreements = -1
		for codon in range(codons):
			if kinetics_codons[codon] > sequence_codons[codon]:
				disagreements += kinetics_codons[codon] - sequence_codons[codon]
				if disagreements >= r:
					break

		# Identify free tRNA candidates
		candidates = 0
		for i in range(trnas_to_codons.shape[1]):
			if trnas_to_codons[codon, i] == 1:
				candidates += free_trnas[i]

		if candidates == 0:
			# If no free tRNAs are available, then the final codon
			# reading event is undone by:
			#	1) Undo the final charging event (returns a charged tRNA
			#	to its free form)
			#	2) Undo the final codon reading event (returns that free
			#	tRNA to its charged form)
			# Net impact: One charging event is undone. One codon
			# reading event is undone. No change to the abundances of
			# free and charged tRNAs.

			# Determine the number of charged tRNAs
			# Assumption: There must be charged tRNAs, since tRNAs exist
			# in either the charged or free form, and there are no free
			# tRNAs.
			for i in range(trnas_to_codons.shape[1]):
				if trnas_to_codons[codon, i] == 1:
					candidates += charged_trnas[i]

			# Randomly select a tRNA that reads this codon
			# Note: prioritizes larger charged tRNA pools
			r = rand() % candidates
			candidates = -1
			for i in range(trnas_to_codons.shape[1]):
				if trnas_to_codons[codon, i] == 1:
					candidates += charged_trnas[i]
					if candidates >= r:
						break

			# Update Kinetics Model
			chargings[i] -= 1
			amino_acids_used[trnas_to_amino_acid_indexes[i]] -= 1
			codons_to_trnas_counter[i, codon] -= 1

		else:
			# If free tRNAs are available, then one free tRNA is
			# returned to its charged form (ie. undo-ing the final codon
			# reading event in the Kinetic Model).

			# Randomly select a tRNA that reads this codon
			# Note: prioritizes larger charged tRNA pools
			r = rand() % candidates
			candidates = -1
			for i in range(trnas_to_codons.shape[1]):
				if trnas_to_codons[codon, i] == 1:
					candidates += free_trnas[i]
					if candidates >= r:
						break

			# Update Kinetics Model
			free_trnas[i] -= 1
			charged_trnas[i] += 1
			codons_to_trnas_counter[i, codon] -= 1

		# Update Kinetics Model
		kinetics_codons[codon] -= 1

	return



@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
@cython.cdivision(True)		# Deactivate 0 division checking
cpdef int get_elongation_rate(
	char [:, :] sequences, # 8-bit signed integer
	int col, # 32-bit signed integer
	double time,
	double target, # 64-bit float
	):

	# Initialize attempted_steps
	cdef int index = 0
	cdef int allocation = 8
	cdef int *attempted_steps = <int *>malloc(
		allocation * sizeof(int))
	if not attempted_steps:
		raise MemoryError()

	# Initialize resulting_rates
	cdef double *resulting_rates = <double *>malloc(
		allocation * sizeof(double))
	if not resulting_rates:
		raise MemoryError()

	# Initialize search parameters
	cdef int lower = 0
	cdef int upper = sequences.shape[1]
	cdef double ribosomes = sequences.shape[0]
	cdef double elongations = 0
	cdef double rate = 0
	cdef double min_difference = 1000 # a large number
	cdef int i = 0
	cdef int j = 0
	cdef bint search_complete = False
	cdef int step_floor = 0
	cdef int step_ceil = 0
	cdef double rate_floor = -1
	cdef double rate_ceil = -1
	cdef int time_int = <int>time
	cdef int result = 0
	cdef int cols = sequences.shape[0]

	# Binary search
	while not search_complete:

		# Count the number of elongations
		elongations = 0
		for i in range(cols):
			for j in range(col):
				if sequences[i, j] != -1:
					elongations += 1

		# Calculate elongation rate
		rate = (elongations
			/ ribosomes
			/ time)

		# Check if more memory is needed
		if index == allocation:

			# Double the memory
			allocation *= 2

			# Reallocate the memory
			attempted_steps = <int*>realloc(
				attempted_steps,
				allocation * sizeof(int))
			resulting_rates = <double*>realloc(
				resulting_rates,
				allocation * sizeof(double))

			# Check that memory has been successfully allocated
			if not attempted_steps:
				raise MemoryError()
			if not resulting_rates:
				raise MemoryError()

		attempted_steps[index] = col
		resulting_rates[index] = rate

		if abs(rate - target) < min_difference:
			min_difference = abs(rate - target)

		if rate > target:
			upper = col
			col = lower + (col - lower) // 2
		elif rate < target:
			lower = col
			col = col + (upper - col) // 2
		else:
			search_complete = True
			break

		index += 1
		for i in range(index):
			if attempted_steps[i] == col:
				search_complete = True
				break

	# Identify the number of columns that resulted in the elongation
	# rate closest to the target value.
	for i in range(index):
		if abs(resulting_rates[i] - target) == min_difference:
			break

	if attempted_steps[i] % 2 == 0:
		result = attempted_steps[i] // time_int

	else:

		step_floor = <int>floor(attempted_steps[i] / time)
		step_ceil = <int>ceil(attempted_steps[i] / time)

		for i in range(index):
			if attempted_steps[i] == step_floor * time_int:
				rate_floor = resulting_rates[i]

			if attempted_steps[i] == step_ceil * time_int:
				rate_ceil = resulting_rates[i]

		if rate_floor == -1:
			# Count the number of elongations
			elongations = 0
			for i in range(cols):
				for j in range(step_floor * time_int):
					if sequences[i, j] != -1:
						elongations += 1

			# Calculate elongation rate
			rate_floor = (elongations
				/ ribosomes
				/ time)

		if rate_ceil == -1:
			# Count the number of elongations
			elongations = 0
			for i in range(cols):
				for j in range(step_ceil * time_int):
					if sequences[i, j] != -1:
						elongations += 1

			# Calculate elongation rate
			rate_ceil = (elongations
				/ ribosomes
				/ time)

		# Identify the number of columns that resulted in the elongation
		# rate closest to the target value.
		if abs(rate_floor - target) < abs(rate_ceil - target):
			result = step_floor
		else:
			result = step_ceil

	free(attempted_steps)
	free(resulting_rates)

	return result


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
@cython.cdivision(True)		# Deactivate 0 division checking
cpdef get_codons_read(
	char [:, :] sequences, # 8 bit signed integer
	long [:] elongations, # 64 bit signed integer
	int size, # 32 bit signed integer
	):

	cdef int i = 0
	cdef int j = 0
	cdef np.ndarray[np.int64_t, ndim=1] out = np.zeros(size, dtype=np.int64)
	cdef np.int64_t[::1] out_view = out

	for i in range(elongations.shape[0]):
		for j in range(elongations[i]):
			out_view[sequences[i, j]] += 1

	return out

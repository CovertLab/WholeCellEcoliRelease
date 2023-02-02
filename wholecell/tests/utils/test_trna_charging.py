"""
Test trna_charging.py

	cd wcEcoli
	pytest wholecell/tests/utils/test_trna_charging.py
"""

from __future__ import absolute_import, division, print_function

from wholecell.utils._trna_charging import (get_initiations,
	get_codon_at, get_candidates_to_C, get_candidates_to_N, select_candidate,
	reconcile_via_ribosome_positions, reconcile_via_trna_pools,
	get_elongation_rate, get_codons_read)

import numpy as np
from numpy.testing import assert_equal

import unittest
from six.moves import range


class Test_trna_charging(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
		pass

	def tearDown(self):
		pass

	def test_get_initiations(self):
		sequence_elongations = np.array([1, 1, 2], dtype=np.int64)
		peptide_lengths = np.array([1, 0, 0], dtype=np.int64)
		protein_indexes = np.array([0, 1, 2], dtype=np.int64)

		n_initiations = get_initiations(
			sequence_elongations,
			peptide_lengths,
			protein_indexes,
			)
		self.assertEqual(n_initiations, 2)

	def test_get_codon_at(self):
		sequences = np.array([
			[0, 1, 2],
			[3, 3, 3],
			], dtype=np.int8)
		elongations = np.array([2, 0], dtype=np.int64)
		ith_ribosome = 0

		# Current codon
		out = get_codon_at(sequences, elongations, ith_ribosome, 0, 0)
		self.assertEqual(out, 1)

		# Codon at +1 position
		out = get_codon_at(sequences, elongations, ith_ribosome, 1, 0)
		self.assertEqual(out, 2)

		# Codon at -1 position
		out = get_codon_at(sequences, elongations, ith_ribosome, -1, 0)
		self.assertEqual(out, 0)

		# Codon beyond right-end (C-terminal) of sequence
		out = get_codon_at(sequences, elongations, ith_ribosome, 2, 0)
		self.assertEqual(out, -1)

		# Codon beyond left-end (N-terminal) of sequence
		out = get_codon_at(sequences, elongations, ith_ribosome, -2, 0)
		self.assertEqual(out, -1)

	def test_get_candidates_to_C(self):
		sequences = np.array([
			[0, 1, 2, 3],
			[0, 2, 1, 3],
			[0, 2, 1, 3],
			], dtype=np.int8)
		elongations = np.array([1, 1, 2], dtype=np.int64)
		candidates = 0
		relative_position = 0
		i = 0

		# Candidates exist in the immediately next (+1) position
		candidates, relative_position = get_candidates_to_C(
			sequences, elongations, 1, candidates, relative_position, i, 0)
		self.assertEqual(candidates, 2)
		self.assertEqual(relative_position, 1)

		# Candidates exist beyond the +1 position
		candidates, relative_position = get_candidates_to_C(
			sequences, elongations, 3, candidates, relative_position, i, 0)
		self.assertEqual(candidates, 1)
		self.assertEqual(relative_position, 2)

		# Candidates do not exist
		candidates, relative_position = get_candidates_to_C(
			sequences, elongations, 0, candidates, relative_position, i, 0)
		self.assertEqual(candidates, 0)
		self.assertEqual(relative_position, 4)

	def test_get_candidates_to_N(self):
		sequences = np.array([
			[0, 1, 2, 3],
			[0, 2, 1, 3],
			[0, 2, 1, 3],
			], dtype=np.int8)
		elongations = np.array([3, 4, 4], dtype=np.int64)
		candidates = 0
		relative_position = 0
		i = 0

		# Candidates exist in the current (0) position
		candidates, relative_position = get_candidates_to_N(
			sequences, elongations, 2, candidates, relative_position, i, 0)
		self.assertEqual(candidates, 1)
		self.assertEqual(relative_position, 0)

		# Candidates exist beyond the 0 position
		candidates, relative_position = get_candidates_to_N(
			sequences, elongations, 1, candidates, relative_position, i, 0)
		self.assertEqual(candidates, 3)
		self.assertEqual(relative_position, -1)

	def test_select_candidate(self):
		sequences = np.array([
			[0, 1, 2, 3],
			[0, 2, 1, 3],
			[0, 2, 1, 3],
			], dtype=np.int8)
		elongations = np.array([1, 1, 2], dtype=np.int64)
		relative_position = 1
		codon_id = 1
		r = 0
		i = 0
		j = 0
		out = select_candidate(
			sequences, elongations, relative_position, codon_id, r, i, j, 0)
		self.assertEqual(out, 0)

		r = 1
		out = select_candidate(
			sequences, elongations, relative_position, codon_id, r, i, j, 0)
		self.assertEqual(out, 2)

	def make_sequence_codons(self, sequences, elongations):
		sequence_codons = np.zeros(sequences.max() + 1, dtype=np.int64)
		for row, cols in enumerate(elongations):
			for col in range(cols):
				sequence_codons[sequences[row, col]] += 1
		return sequence_codons

	def make_codons_to_trnas_counter(self, free_trnas, sequence_codons):
		n_trnas = free_trnas.shape[0]
		n_codons = sequence_codons.shape[0]
		codons_to_trnas_counter = np.zeros((n_trnas, n_codons), dtype=np.int64)
		return codons_to_trnas_counter

	def test_reconcile_equal(self):
		# No disagreements to reconcile
		kinetics_codons = np.array([0, 1, 0], dtype=np.int64)
		elongations = np.array([1, 0, 0], dtype=np.int64)
		sequences = np.array([
			[1, 0, 1, 2],
			[1, 2, 1, 2],
			[1, 1, 2, 0],
			], dtype=np.int8)
		sequence_codons = self.make_sequence_codons(sequences, elongations)

		reconcile_via_ribosome_positions(
			sequence_codons,
			elongations,
			kinetics_codons,
			sequences,
			np.byte(4),
			)

		assert_equal(sequence_codons, np.array([0, 1, 0], dtype=np.int64))
		assert_equal(elongations, np.array([1, 0, 0], dtype=np.int64))
		assert_equal(kinetics_codons, np.array([0, 1, 0], dtype=np.int64))

	def test_reconcile_forward(self):
		# Need to take +1 step to reconcile
		kinetics_codons = np.array([0, 2, 1], dtype=np.int64)
		elongations = np.array([1, 1, 0], dtype=np.int64)
		sequences = np.array([
			[1, 0, 1, 2],
			[1, 2, 1, 2],
			[1, 1, 2, 0],
			], dtype=np.int8)
		sequence_codons = self.make_sequence_codons(sequences, elongations)

		reconcile_via_ribosome_positions(
			sequence_codons,
			elongations,
			kinetics_codons,
			sequences,
			np.byte(4),
			)

		assert_equal(sequence_codons, np.array([0, 2, 1], dtype=np.int64))
		assert_equal(elongations, np.array([1, 2, 0], dtype=np.int64))
		assert_equal(kinetics_codons, np.array([0, 2, 1], dtype=np.int64))

	def test_reconcile_backward(self):
		# Need to take -1 step to reconcile
		kinetics_codons = np.array([1, 6, 1], dtype=np.int64)
		elongations = np.array([3, 3, 3], dtype=np.int64)
		sequences = np.array([
			[1, 0, 1, 2],
			[1, 2, 1, 2],
			[1, 1, 2, 0],
			], dtype=np.int8)
		sequence_codons = self.make_sequence_codons(sequences, elongations)

		reconcile_via_ribosome_positions(
			sequence_codons,
			elongations,
			kinetics_codons,
			sequences,
			np.byte(4),
			)

		assert_equal(sequence_codons, np.array([1, 6, 1], dtype=np.int64))
		assert_equal(elongations, np.array([3, 3, 2], dtype=np.int64))
		assert_equal(kinetics_codons, np.array([1, 6, 1], dtype=np.int64))

	def test_reconcile_use_free_trna(self):
		# Reconciling using only ribosome positions cannot resolve
		kinetics_codons = np.array([0, 0, 1, 0], dtype=np.int64)
		elongations = np.array([0, 0, 0], dtype=np.int64)
		sequences = np.array([
			[1, 0, 1, 2],
			[3, 2, 1, 2],
			[1, 1, 2, 0],
			], dtype=np.int8)
		sequence_codons = self.make_sequence_codons(sequences, elongations)
		
		reconcile_via_ribosome_positions(
			sequence_codons,
			elongations,
			kinetics_codons,
			sequences,
			np.byte(4),
			)

		assert_equal(sequence_codons, np.array([0, 0, 0, 0], dtype=np.int64))
		assert_equal(elongations, np.array([0, 0, 0], dtype=np.int64))
		assert_equal(kinetics_codons, np.array([0, 0, 1, 0], dtype=np.int64))

		# Need to use 1 free trna to reconcile
		free_trnas = np.array([2, 0], dtype=np.int64)
		charged_trnas = np.array([0, 0], dtype=np.int64)
		chargings = np.array([1, 0], dtype=np.int64)
		amino_acids_used = np.array([1], dtype=np.int64)
		trnas_to_amino_acid_indexes = np.array([0, 0], dtype=np.int8)
		trnas_to_codons = np.array([
			[1, 0],
			[1, 1],
			[1, 0],
			[1, 1],
			], dtype=np.int8)
		codons_to_trnas_counter = self.make_codons_to_trnas_counter(free_trnas, sequence_codons)

		# tRNA 0 interacted with codon 2, once
		codons_to_trnas_counter[0, 2] = 1

		reconcile_via_trna_pools(
			sequence_codons,
			kinetics_codons,
			free_trnas,
			charged_trnas,
			chargings,
			amino_acids_used,
			codons_to_trnas_counter,
			trnas_to_codons,
			trnas_to_amino_acid_indexes,
			)
		
		assert_equal(free_trnas, np.array([1, 0], dtype=np.int64))
		assert_equal(charged_trnas, np.array([1, 0], dtype=np.int64))
		assert_equal(chargings, np.array([1, 0], dtype=np.int64))
		assert_equal(amino_acids_used, np.array([1], dtype=np.int64))
		assert_equal(
			codons_to_trnas_counter,
			np.array([[0, 0, 0, 0], [0, 0, 0, 0]], dtype=np.int64))

	def test_reconcile_forward_undo_charging(self):
		# No candidates
		kinetics_codons = np.array([0, 0, 1, 0], dtype=np.int64)
		elongations = np.array([0, 0, 0], dtype=np.int64)
		sequences = np.array([
			[1, 0, 1, 2],
			[3, 2, 1, 2],
			[1, 1, 2, 0],
			], dtype=np.int8)
		sequence_codons = self.make_sequence_codons(sequences, elongations)
		
		reconcile_via_ribosome_positions(
			sequence_codons,
			elongations,
			kinetics_codons,
			sequences,
			np.byte(4),
			)

		assert_equal(sequence_codons, np.array([0, 0, 0, 0], dtype=np.int64))
		assert_equal(elongations, np.array([0, 0, 0], dtype=np.int64))
		assert_equal(kinetics_codons, np.array([0, 0, 1, 0], dtype=np.int64))

		# Need to undo 1 charging event, then use 1 free trna to reconcile
		free_trnas = np.array([0, 0], dtype=np.int64)
		charged_trnas = np.array([2, 0], dtype=np.int64)
		chargings = np.array([1, 0], dtype=np.int64)
		amino_acids_used = np.array([1], dtype=np.int64)
		trnas_to_amino_acid_indexes = np.array([0, 0], dtype=np.int8)
		trnas_to_codons = np.array([
			[1, 0],
			[1, 1],
			[1, 0],
			[1, 1],
			], dtype=np.int8)
		codons_to_trnas_counter = self.make_codons_to_trnas_counter(free_trnas, sequence_codons)

		# tRNA 0 interacted with codon 2, once
		codons_to_trnas_counter[0, 2] = 1

		reconcile_via_trna_pools(
			sequence_codons,
			kinetics_codons,
			free_trnas,
			charged_trnas,
			chargings,
			amino_acids_used,
			codons_to_trnas_counter,
			trnas_to_codons,
			trnas_to_amino_acid_indexes,
			)

		assert_equal(free_trnas, np.array([0, 0], dtype=np.int64))
		assert_equal(charged_trnas, np.array([2, 0], dtype=np.int64))
		assert_equal(chargings, np.array([0, 0], dtype=np.int64))
		assert_equal(amino_acids_used, np.array([0], dtype=np.int64))
		assert_equal(
			codons_to_trnas_counter,
			np.array([[0, 0, 0, 0], [0, 0, 0, 0]], dtype=np.int64))

	def test_reconcile_backward_beyond(self):
		# Need to take more than -1 step to reconcile
		kinetics_codons = np.array([3, 0, 3, 0], dtype=np.int64)
		elongations = np.array([3, 3, 3], dtype=np.int64)
		sequences = np.array([
			[0, 0, 2, 2],
			[0, 2, 2, 2],
			[1, 3, 3, 0],
			], dtype=np.int8)
		sequence_codons = self.make_sequence_codons(sequences, elongations)
		
		reconcile_via_ribosome_positions(
			sequence_codons,
			elongations,
			kinetics_codons,
			sequences,
			np.byte(4),
			)

		assert_equal(sequence_codons, np.array([3, 0, 3, 0], dtype=np.int64))
		assert_equal(kinetics_codons, np.array([3, 0, 3, 0], dtype=np.int64))
		assert_equal(elongations, np.array([3, 3, 0], dtype=np.int64))

	def test_reconcile_attempts_threshold(self):
		kinetics_codons = np.array([10, 20], dtype=np.int64)
		elongations = 3 * np.ones(10, dtype=np.int64)
		sequences = np.array([
			[0, 1, 0, 1, 0],
			[1, 0, 1, 0, 1],
			[0, 1, 0, 1, 0],
			[1, 0, 1, 0, 1],
			[0, 1, 0, 1, 0],
			[1, 0, 1, 0, 1],
			[0, 1, 0, 1, 0],
			[1, 0, 1, 0, 1],
			[0, 1, 0, 1, 0],
			[1, 0, 1, 0, 1],
			], dtype=np.int8)
		sequence_codons = self.make_sequence_codons(sequences, elongations)
		
		reconcile_via_ribosome_positions(
			sequence_codons,
			elongations,
			kinetics_codons,
			sequences,
			np.byte(4),
			)

		assert_equal(sequence_codons, np.array([10, 15], dtype=np.int64))
		assert_equal(kinetics_codons, np.array([10, 20], dtype=np.int64))
		self.assertEqual(elongations.sum(), 25)

	def test_get_elongation_rate(self):
		sequences = np.array([
			[0, 1, 1, -1, -1],
			[1, 0, 0, 1, 1],
			],
			dtype=np.int8)
		previous_rate = 3
		time_step = 1
		target = 4

		rate = get_elongation_rate(
			sequences,
			previous_rate,
			time_step,
			target
			)
		assert_equal(rate, 4)

	def test_get_codons_read(self):
		sequences = np.array([
			[0, 1, 1, -1, -1],
			[1, 0, 0, 1, 1],
			],
			dtype=np.int8)
		elongations = np.array([2, 4], dtype=np.int64)
		size = 2

		n_codons_read = get_codons_read(
			sequences,
			elongations,
			size
			)
		assert_equal(n_codons_read, np.array([3, 3], dtype=np.int64))


if __name__ == '__main__':
	unittest.main()


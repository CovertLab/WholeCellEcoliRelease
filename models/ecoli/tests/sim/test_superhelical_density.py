"""
test_superhelical_density.py

@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 4/28/2020
"""
from __future__ import absolute_import, division, print_function

import unittest

import numpy as np
import numpy.testing as npt

from models.ecoli.processes.chromosome_structure import ChromosomeStructure


class TestSuperhelicalDensity(unittest.TestCase):
	"""
	Tests the method self._compute_new_segment_attributes() defined in the
	process ChromosomeStructure.
	"""
	def setUp(self):
		self.terC_index = -1
		self.min_coordinates = -12
		self.max_coordinates = 12

		self.cs = ChromosomeStructure()
		self.cs.min_coordinates = self.min_coordinates
		self.cs.max_coordinates = self.max_coordinates
		self.cs.terC_index = self.terC_index

		self.compute_new_segment_attributes = self.cs._compute_new_segment_attributes

	def tearDown(self):
		pass


	def test_no_change(self):
		"""
		Test all attributes are preserved when there are no changes.
		"""
		boundary_molecule_indexes = np.array([[0, 1], [1, 2]])
		boundary_coordinates = np.array([[-8, -2], [-2, 4]])
		linking_numbers = np.array([3, 2])

		new_molecule_indexes = np.array([0, 1, 2])
		new_molecule_coordinates = np.array([-8, -2, 4])
		domain_spans_oriC = True
		domain_spans_terC = False

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes, new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], linking_numbers)
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			boundary_molecule_indexes)
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			boundary_coordinates)

		# Shuffled inputs
		boundary_molecule_indexes_shuffled = np.array([[1, 2], [0, 1]])
		boundary_coordinates_shuffled = np.array([[-2, 4], [-8, -2]])
		linking_numbers_shuffled = np.array([2, 3])

		new_molecule_indexes_shuffled = np.array([2, 1, 0])
		new_molecule_coordinates_shuffled = np.array([4, -2, -8])

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes_shuffled, boundary_coordinates_shuffled,
			linking_numbers_shuffled, new_molecule_indexes_shuffled,
			new_molecule_coordinates_shuffled, domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], linking_numbers)
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			boundary_molecule_indexes)
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			boundary_coordinates)


	def test_constant_linking_number(self):
		"""
		Test linking numbers are preserved when no new molecules are added or
		removed.
		"""
		boundary_molecule_indexes = np.array([[0, 1], [1, 2]])
		boundary_coordinates = np.array([[-8, -2], [-2, 4]])
		linking_numbers = np.array([3, 2])

		new_molecule_indexes = np.array([0, 1, 2])
		new_molecule_coordinates = np.array([-6, -4, 8])
		domain_spans_oriC = True
		domain_spans_terC = False

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes, new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], linking_numbers)
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			boundary_molecule_indexes)
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[-6, -4], [-4, 8]]))


	def test_molecules_added(self):
		"""
		Test linking numbers are split correctly when new molecules are added.
		"""
		boundary_molecule_indexes = np.array([[0, 1], [1, 2]])
		boundary_coordinates = np.array([[-8, -2], [-2, 4]])
		linking_numbers = np.array([3, 3])

		new_molecule_indexes = np.array([0, 1, 2, 3, 4])
		new_molecule_coordinates = np.array([-8, -2, 4, -4, 2])
		domain_spans_oriC = True
		domain_spans_terC = False

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes, new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], np.array([2, 1, 2, 1]))
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			np.array([[0, 3], [3, 1], [1, 4], [4, 2]]))
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[-8, -4], [-4, -2], [-2, 2], [2, 4]]))


	def test_molecules_removed(self):
		"""
		Test linking numbers are added up correctly when molecules are removed.
		"""
		boundary_molecule_indexes = np.array([[0, 1], [1, 2]])
		boundary_coordinates = np.array([[-8, -2], [-2, 4]])
		linking_numbers = np.array([3, 2])

		new_molecule_indexes = np.array([0, 2])
		new_molecule_coordinates = np.array([-8, 4])
		domain_spans_oriC = True
		domain_spans_terC = False

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes, new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], np.array([5]))
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			np.array([[0, 2]]))
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[-8, 4]]))


	def test_molecules_added_and_removed(self):
		"""
		Test linking numbers are calculated correctly when molecules are
		simultaneously added and removed.
		"""
		boundary_molecule_indexes = np.array([[0, 1], [1, 2]])
		boundary_coordinates = np.array([[-8, -2], [-2, 4]])
		linking_numbers = np.array([4, 2])

		new_molecule_indexes = np.array([0, 2, 3])
		new_molecule_coordinates = np.array([-8, 4, -4])
		domain_spans_oriC = True
		domain_spans_terC = False

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes, new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], np.array([2, 4]))
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			np.array([[0, 3], [3, 2]]))
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[-8, -4], [-4, 4]]))


	def test_terC_segment(self):
		"""
		Test that segments spanning terC are handled properly.
		"""
		boundary_molecule_indexes = np.array([[self.terC_index, 0], [0, 1], [2, self.terC_index]])
		boundary_coordinates = np.array([[self.min_coordinates, -8], [-8, -4], [6, self.max_coordinates]])
		linking_numbers = np.array([2, 1, 3])

		# No molecules added or removed
		new_molecule_indexes = np.array([0, 1, 2])
		new_molecule_coordinates = np.array([-6, -2, 6])
		domain_spans_oriC = False
		domain_spans_terC = True

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes, new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], np.array([2.5, 1, 2.5]))
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			boundary_molecule_indexes)
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[self.min_coordinates, -6], [-6, -2], [6, self.max_coordinates]]))

		# Molecules removed
		new_molecule_indexes = np.array([1, 2])
		new_molecule_coordinates = np.array([-6, 6])
		domain_spans_oriC = False
		domain_spans_terC = True

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes, new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], np.array([3, 3]))
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			np.array([[self.terC_index, 1], [2, self.terC_index]]))
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[self.min_coordinates, -6], [6, self.max_coordinates]]))

		# Molecules added
		new_molecule_indexes = np.array([0, 1, 2, 3])
		new_molecule_coordinates = np.array([-6, -2, 6, 9])
		domain_spans_oriC = False
		domain_spans_terC = True

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes, new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], np.array([2.5, 1, 1.25, 1.25]))
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			np.array([[self.terC_index, 0], [0, 1], [2, 3], [3, self.terC_index]]))
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[self.min_coordinates, -6], [-6, -2], [6, 9], [9, self.max_coordinates]]))


	def test_empty_chromosome(self):
		"""
		Test handling of empty chromosomes with no bound molecules.
		"""
		# Empty to empty
		boundary_molecule_indexes = np.array([[self.terC_index, self.terC_index]])
		boundary_coordinates = np.array([[self.min_coordinates, self.max_coordinates]])
		linking_numbers = np.array([24])

		new_molecule_indexes = np.array([])
		new_molecule_coordinates = np.array([])
		domain_spans_oriC = True
		domain_spans_terC = True

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes,
			new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], linking_numbers)
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			np.array([[self.terC_index, self.terC_index]]))
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[self.min_coordinates, self.max_coordinates]]))

		# Empty to nonempty
		new_molecule_indexes = np.array([0, 1])
		new_molecule_coordinates = np.array([-2, 4])
		domain_spans_oriC = True
		domain_spans_terC = True

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes,
			new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], np.array([10, 6, 8]))
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			np.array([[self.terC_index, 0], [0, 1], [1, self.terC_index]]))
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[self.min_coordinates, -2], [-2, 4], [4, self.max_coordinates]]))

		# Nonempty to empty
		boundary_molecule_indexes = np.array([[self.terC_index, 0], [0, 1], [1, self.terC_index]])
		boundary_coordinates = np.array([[self.min_coordinates, -2], [-2, 4], [4, self.max_coordinates]])
		linking_numbers = np.array([10, 6, 8])

		new_molecule_indexes = np.array([])
		new_molecule_coordinates = np.array([])
		domain_spans_oriC = True
		domain_spans_terC = True

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes,
			new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], np.array([24]))
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			np.array([[self.terC_index, self.terC_index]]))
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[self.min_coordinates, self.max_coordinates]]))


	def test_initiation(self):
		"""
		Test handling of case where the segment spanning oriC becomes
		disconnected (new round of replication initiated). Note that the two
		replisomes are always initiated on the origin (coordinate zero).
		"""
		boundary_molecule_indexes = np.array([[0, 1], [1, 2], [2, 3]])
		boundary_coordinates = np.array([[-10, -4], [-4, 4], [4, 8]])
		linking_numbers = np.array([1.5, 2, 1])

		new_molecule_indexes = np.array([0, 1, 2, 3, 4, 5])
		new_molecule_coordinates = np.array([-10, -4, 4, 8, 0, 0])
		domain_spans_oriC = False
		domain_spans_terC = False

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes,
			new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], np.array([1.5, 1, 1, 1]))
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			np.array([[0, 1], [1, 4], [5, 2], [2, 3]]))
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[-10, -4], [-4, 0], [0, 4], [4, 8]]))


	def test_termination(self):
		"""
		Test handling of case where the given chromosomal domain connects at
		terC (round of replication completed).
		"""
		boundary_molecule_indexes = np.array([[0, 1], [2, 3]])
		boundary_coordinates = np.array([[-10, -4], [6, 10]])
		linking_numbers = np.array([3, 4])

		new_molecule_indexes = np.array([1, 2])
		new_molecule_coordinates = np.array([-5, 5])
		domain_spans_oriC = False
		domain_spans_terC = True

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes,
			new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], np.array([3.5, 3.5]))
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			np.array([[self.terC_index, 1], [2, self.terC_index]]))
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[self.min_coordinates, -5], [5, self.max_coordinates]]))


	def test_simultaneous_initiation_termination(self):
		"""
		Test handling of case where the initiation of a new round of
		replication and the termination of the current round of replication
		occurs simultaneously.
		"""
		boundary_molecule_indexes = np.array([[0, 1], [1, 2], [2, 3]])
		boundary_coordinates = np.array([[-10, -4], [-4, 6], [6, 10]])
		linking_numbers = np.array([3, 2, 4])

		new_molecule_indexes = np.array([1, 2, 4, 5])
		new_molecule_coordinates = np.array([-5, 5, 0, 0])
		domain_spans_oriC = False
		domain_spans_terC = True

		new_segment_attrs = self.compute_new_segment_attributes(
			boundary_molecule_indexes, boundary_coordinates,
			linking_numbers, new_molecule_indexes,
			new_molecule_coordinates,
			domain_spans_oriC, domain_spans_terC
			)

		npt.assert_equal(
			new_segment_attrs['linking_numbers'], np.array([3.5, 1, 1, 3.5]))
		npt.assert_equal(
			new_segment_attrs['boundary_molecule_indexes'],
			np.array([[self.terC_index, 1], [1, 4], [5, 2], [2, self.terC_index]]))
		npt.assert_equal(
			new_segment_attrs['boundary_coordinates'],
			np.array([[self.min_coordinates, -5], [-5, 0], [0, 5], [5, self.max_coordinates]]))

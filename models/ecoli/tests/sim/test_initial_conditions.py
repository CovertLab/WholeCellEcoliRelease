#!/usr/bin/env python
"""
Test_InitialConditions.py
Tests the initial conditions code of the model.
"""

from __future__ import absolute_import, division, print_function

import unittest

import numpy as np
import six

from wholecell.utils import units

from  models.ecoli.sim.initial_conditions import determine_chromosome_state


N_MAX_REPLISOMES = 1000  # Chosen to be big enough to have no effect


class Test_InitialConditions(unittest.TestCase):
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


	def test_num_forks(self):

		# When (C + D) / tau is less than one, no replication will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 90. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		expected_num_forks = 0

		self.assertEqual(expected_num_forks, len(replisome_state["coordinates"]))
		self.assertEqual(expected_num_forks, len(replisome_state["right_replichore"]))
		self.assertEqual(expected_num_forks, len(replisome_state["domain_index"]))

		# When (C + D) / tau is between one and two, one replication generation will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 60. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		expected_num_forks = 2

		self.assertEqual(expected_num_forks, len(replisome_state["coordinates"]))
		self.assertEqual(expected_num_forks, len(replisome_state["right_replichore"]))
		self.assertEqual(expected_num_forks, len(replisome_state["domain_index"]))

		# When (C + D) / tau is between two and three, two replication generations will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 30. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		# Two generations means one event from first generation and two events
		# from the second - total of 3 events, 2 forks per event = 6 replisomes
		expected_num_forks = 6

		self.assertEqual(expected_num_forks, len(replisome_state["coordinates"]))
		self.assertEqual(expected_num_forks, len(replisome_state["right_replichore"]))
		self.assertEqual(expected_num_forks, len(replisome_state["domain_index"]))

		# When (C + D) / tau is between three and four, three replication generations will have started
		C = 50. * units.min
		D = 19. * units.min
		tau = 20. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		# Three generations means one event from first generation, two from the
		# second, and 4 from the third - total of 7 events,
		# 2 forks per event = 14 replisomes
		expected_num_forks = 14

		self.assertEqual(expected_num_forks, len(replisome_state["coordinates"]))
		self.assertEqual(expected_num_forks, len(replisome_state["right_replichore"]))
		self.assertEqual(expected_num_forks, len(replisome_state["domain_index"]))

		# When (C + D) / tau is between four and five, four replication
		# generations will have started
		C = 70. * units.min
		D = 20. * units.min
		tau = 21. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		# Four generations means one event from first generation, two from the
		# second, 4 from the third, and 8 from the fourth - total of 15 events,
		# 2 forks per event = 30 replisomes
		expected_num_forks = 30

		self.assertEqual(expected_num_forks, len(replisome_state["coordinates"]))
		self.assertEqual(expected_num_forks, len(replisome_state["right_replichore"]))
		self.assertEqual(expected_num_forks, len(replisome_state["domain_index"]))


	def test_fork_coordinates(self):

		C = 70. * units.min
		D = 25. * units.min
		tau = 30. * units.min
		replichore_length = 2319838 * units.nt

		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)
		fork_coordinates = replisome_state["coordinates"]

		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.strip_empty_units(ratio)
			expected_coordinates = np.floor(
				ratio * (replichore_length.asNumber()))
			self.assertEqual(expected_coordinates,
				fork_coordinates[2 ** n - 2])
			self.assertEqual(expected_coordinates,
				-fork_coordinates[2 ** n - 1])
			n += 1

		# Length of DNA to be replicated is 1
		C = 70. * units.min
		D = 25. * units.min
		tau = 30. * units.min
		replichore_length = 1 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)
		fork_coordinates = replisome_state["coordinates"]

		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.strip_empty_units(ratio)
			expected_coordinates = np.floor(
				ratio * (replichore_length.asNumber()))
			self.assertEqual(expected_coordinates,
				fork_coordinates[2 ** n - 2])
			self.assertEqual(expected_coordinates,
				-fork_coordinates[2 ** n - 1])
			n += 1

		# Length of DNA to be replicated is 0
		C = 70. * units.min
		D = 25. * units.min
		tau = 30. * units.min
		replichore_length = 0 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)
		fork_coordinates = replisome_state["coordinates"]

		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.strip_empty_units(ratio)
			expected_coordinates = np.floor(
				ratio * (replichore_length.asNumber()))
			self.assertEqual(expected_coordinates,
				fork_coordinates[2 ** n - 2])
			self.assertEqual(expected_coordinates,
				-fork_coordinates[2 ** n - 1])
			n += 1

		# Length of DNA to be replicated is not a whole number
		C = 70. * units.min
		D = 25. * units.min
		tau = 30. * units.min
		replichore_length = 2319838.3 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)
		fork_coordinates = replisome_state["coordinates"]

		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.strip_empty_units(ratio)
			expected_coordinates = np.floor(
				ratio * (replichore_length.asNumber()))
			self.assertEqual(expected_coordinates,
				fork_coordinates[2 ** n - 2])
			self.assertEqual(expected_coordinates,
				-fork_coordinates[2 ** n - 1])
			n += 1


	def test_maximum_domain_index(self):

		# When (C + D) / tau is less than one, no replication will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 90. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		self.assertEqual(0, domain_state["domain_index"].max())


		# When (C + D) / tau is between one and two, one replication generation will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 60. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		self.assertEqual(2, domain_state["domain_index"].max())


		# When (C + D) / tau is between two and three, two replication generations will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 30. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		# There should be equal numbers of zeros and ones in the chromosomeIndex array (excepting the first four)
		self.assertEqual(6, domain_state["domain_index"].max())


		# When (C + D) / tau is between three and four, three replication generations will have started
		C = 50. * units.min
		D = 19. * units.min
		tau = 20. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		# There should be equal numbers of zeros and ones in the chromosomeIndex array (excepting the first four)
		self.assertEqual(14, domain_state["domain_index"].max())


		# When (C + D) / tau is between four and five, four replication generations will have started
		C = 70. * units.min
		D = 20. * units.min
		tau = 21. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		# There should be equal numbers of zeros and ones in the chromosomeIndex array (excepting the first four)
		self.assertEqual(30, domain_state["domain_index"].max())


	def test_determine_chromosome_state_inputs(self):

		# The D period must be shorter than tau
		with six.assertRaisesRegex(self,
			AssertionError,
			"^The D period must be shorter than the doubling time tau.$"
		):
			C = 50. * units.min
			D = 20. * units.min
			tau = 19. * units.min
			replichore_length = 2319838 * units.nt
			oric_state, replisome_state, domain_state = determine_chromosome_state(
				C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		# No inputs can have negative values
		with six.assertRaisesRegex(self,
			AssertionError,
			"^C value can't be negative.$"
		):
			C = -50. * units.min
			D = 20. * units.min
			tau = 60. * units.min
			replichore_length = 2319838 * units.nt
			oric_state, replisome_state, domain_state = determine_chromosome_state(
				C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		with six.assertRaisesRegex(self,
			AssertionError,
			"^D value can't be negative.$"
		):
			C = 50. * units.min
			D = -20. * units.min
			tau = 60. * units.min
			replichore_length = 2319838 * units.nt
			oric_state, replisome_state, domain_state = determine_chromosome_state(
				C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		with six.assertRaisesRegex(self,
			AssertionError,
			"^tau value can't be negative.$"
		):
			C = 50. * units.min
			D = 20. * units.min
			tau = -60. * units.min
			replichore_length = 2319838 * units.nt
			oric_state, replisome_state, domain_state = determine_chromosome_state(
				C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		with six.assertRaisesRegex(self,
			AssertionError,
			"^replichore_length value can't be negative.$"
		):
			C = 50. * units.min
			D = 20. * units.min
			tau = 60. * units.min
			replichore_length = -2319838 * units.nt
			oric_state, replisome_state, domain_state = determine_chromosome_state(
				C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)


	def test_num_oriCs(self):

		# When (C + D) / tau is less than one, no replication will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 100. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		self.assertEqual(1, len(oric_state["domain_index"]))

		# When (C + D) / tau is less than one, no replication will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 90. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		self.assertEqual(1, len(oric_state["domain_index"]))

		# When (C + D) / tau is between one and two, one replication generation will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 60. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		self.assertEqual(2, len(oric_state["domain_index"]))

		# When (C + D) / tau is between two and three, two replication generations will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 30. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			C, D, tau, replichore_length, N_MAX_REPLISOMES, -1)

		# Two generations means one event from first generation and two events
		# from the second - total of 3 events, plus one oriC for the original
		# chromosome = 4 total oriC
		self.assertEqual(4, len(oric_state["domain_index"]))


if __name__ == '__main__':
	unittest.main()

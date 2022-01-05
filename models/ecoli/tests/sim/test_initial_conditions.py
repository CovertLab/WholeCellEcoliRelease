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

from reconstruction.ecoli.initialization import determine_chromosome_state


N_MAX_REPLISOMES = 1000  # Chosen to be big enough to have no effect
REPLICATION_RATE = 967


def cell_mass_from_tau(tau):
	"""Approximate initial cell mass based on fit from data calculated in
	growth_rate_dependent_parameters.py"""
	return 2.5 / (7e-5 * tau.asNumber(units.min) - 6e-4) * units.fg

CRITICAL_MASS = cell_mass_from_tau(60 * units.min)


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

		# When cell mass / critical mass is less than one, no replication will have started
		tau = 90. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		expected_num_forks = 0

		self.assertEqual(expected_num_forks, len(replisome_state["coordinates"]))
		self.assertEqual(expected_num_forks, len(replisome_state["right_replichore"]))
		self.assertEqual(expected_num_forks, len(replisome_state["domain_index"]))

		# When cell mass / critical mass is between one and two, one replication generation will have started
		tau = 50. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		expected_num_forks = 2

		self.assertEqual(expected_num_forks, len(replisome_state["coordinates"]))
		self.assertEqual(expected_num_forks, len(replisome_state["right_replichore"]))
		self.assertEqual(expected_num_forks, len(replisome_state["domain_index"]))

		# When cell mass / critical mass is between two and four, two replication generations will have started
		tau = 30. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		# Two generations means one event from first generation and two events
		# from the second - total of 3 events, 2 forks per event = 6 replisomes
		expected_num_forks = 6

		self.assertEqual(expected_num_forks, len(replisome_state["coordinates"]))
		self.assertEqual(expected_num_forks, len(replisome_state["right_replichore"]))
		self.assertEqual(expected_num_forks, len(replisome_state["domain_index"]))

		# When cell mass / critical mass is between four and eight, three replication generations will have started
		tau = 20. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		# Three generations means one event from first generation, two from the
		# second, and 4 from the third - total of 7 events,
		# 2 forks per event = 14 replisomes
		expected_num_forks = 14

		self.assertEqual(expected_num_forks, len(replisome_state["coordinates"]))
		self.assertEqual(expected_num_forks, len(replisome_state["right_replichore"]))
		self.assertEqual(expected_num_forks, len(replisome_state["domain_index"]))

		# When cell mass / critical mass is between eight and 16, four replication
		# generations will have started
		tau = 12. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		# Four generations means one event from first generation, two from the
		# second, 4 from the third, and 8 from the fourth - total of 15 events,
		# 2 forks per event = 30 replisomes
		expected_num_forks = 30

		self.assertEqual(expected_num_forks, len(replisome_state["coordinates"]))
		self.assertEqual(expected_num_forks, len(replisome_state["right_replichore"]))
		self.assertEqual(expected_num_forks, len(replisome_state["domain_index"]))


	def test_fork_coordinates(self):

		tau = 30. * units.min
		replichore_length = 2319838 * units.nt

		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)
		fork_coordinates = replisome_state["coordinates"]
		expected_coordinates = np.array([2198437, -2198437, 457837, -457837, 457837, -457837])
		np.testing.assert_array_equal(fork_coordinates, expected_coordinates)

		# Length of DNA to be replicated is 1
		tau = 30. * units.min
		replichore_length = 1 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)
		fork_coordinates = replisome_state["coordinates"]
		expected_coordinates = np.array([0, 0, 0, 0, 0, 0])
		np.testing.assert_array_equal(fork_coordinates, expected_coordinates)

		# Length of DNA to be replicated is not a whole number
		tau = 30. * units.min
		replichore_length = 2319838.3 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)
		fork_coordinates = replisome_state["coordinates"]
		expected_coordinates = np.array([2198437, -2198437, 457837, -457837, 457837, -457837])
		np.testing.assert_array_equal(fork_coordinates, expected_coordinates)


	def test_maximum_domain_index(self):

		# When cell mass / critical mass is less than one, no replication will have started
		tau = 90. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		self.assertEqual(0, domain_state["domain_index"].max())


		# When cell mass / critical mass is between one and two, one replication generation will have started
		tau = 50. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		self.assertEqual(2, domain_state["domain_index"].max())


		# When cell mass / critical mass is between two and four, two replication generations will have started
		tau = 30. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		# There should be equal numbers of zeros and ones in the chromosomeIndex array (excepting the first four)
		self.assertEqual(6, domain_state["domain_index"].max())


		# When cell mass / critical mass is between four and eight, three replication generations will have started
		tau = 20. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		# There should be equal numbers of zeros and ones in the chromosomeIndex array (excepting the first four)
		self.assertEqual(14, domain_state["domain_index"].max())


		# When cell mass / critical mass is between eight and 16, four replication generations will have started
		tau = 12. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		# There should be equal numbers of zeros and ones in the chromosomeIndex array (excepting the first four)
		self.assertEqual(30, domain_state["domain_index"].max())


	def test_determine_chromosome_state_inputs(self):

		with six.assertRaisesRegex(self,
			AssertionError,
			"^tau value can't be negative.$"
		):
			tau = -60. * units.min
			replichore_length = 2319838 * units.nt
			oric_state, replisome_state, domain_state = determine_chromosome_state(
				tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		with six.assertRaisesRegex(self,
			AssertionError,
			"^replichore_length must be positive.$"
		):
			tau = 60. * units.min
			replichore_length = -2319838 * units.nt
			oric_state, replisome_state, domain_state = determine_chromosome_state(
				tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)


	def test_num_oriCs(self):

		# When cell mass / critical mass is less than one, no replication will have started
		tau = 100. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		self.assertEqual(1, len(oric_state["domain_index"]))

		# When cell mass / critical mass is less than one, no replication will have started
		tau = 90. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		self.assertEqual(1, len(oric_state["domain_index"]))

		# When cell mass / critical mass is between one and two, one replication generation will have started
		tau = 50. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		self.assertEqual(2, len(oric_state["domain_index"]))

		# When cell mass / critical mass is between two and four, two replication generations will have started
		tau = 30. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(
			tau, replichore_length, N_MAX_REPLISOMES, -1, cell_mass_from_tau(tau), CRITICAL_MASS, REPLICATION_RATE)

		# Two generations means one event from first generation and two events
		# from the second - total of 3 events, plus one oriC for the original
		# chromosome = 4 total oriC
		self.assertEqual(4, len(oric_state["domain_index"]))


if __name__ == '__main__':
	unittest.main()

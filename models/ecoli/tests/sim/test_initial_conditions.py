"""
Test_InitialConditions.py
Tests the initial conditions code of the model.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/23/2015
"""

from __future__ import division

import unittest
import warnings
import nose.plugins.attrib as noseAttrib

import numpy as np
import cPickle
import os

from wholecell.utils import units

from  models.ecoli.sim.initial_conditions import determine_chromosome_state


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


	@noseAttrib.attr('replicationTest')
	def test_num_forks(self):

		# When (C + D) / tau is less than one, no replication will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 90. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		self.assertTrue(len(replisome_state["coordinates"]) == len(replisome_state["right_replichore"]) == len(replisome_state["domain_index"]) == 0)

		# When (C + D) / tau is between one and two, one replication generation will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 60. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		self.assertTrue(len(replisome_state["coordinates"]) == len(replisome_state["right_replichore"]) == len(replisome_state["domain_index"]) == 2)

		# When (C + D) / tau is between two and three, two replication generations will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 30. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		# Two generations means one event from first generation and two events
		# from the second - total of 3 events, 2 forks per event = 6 replisomes
		self.assertTrue(len(replisome_state["coordinates"]) == len(replisome_state["right_replichore"]) == len(replisome_state["domain_index"]) == 6)

		# When (C + D) / tau is between three and four, three replication generations will have started
		C = 50. * units.min
		D = 19. * units.min
		tau = 20. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		# Three generations means one event from first generation, two from the
		# second, and 4 from the third - total of 7 events, 2 forks per event = 14 replisomes
		self.assertTrue(len(replisome_state["coordinates"]) == len(replisome_state["right_replichore"]) == len(replisome_state["domain_index"]) == 14)

		# When (C + D) / tau is between four and five, four replication generations will have started
		C = 70. * units.min
		D = 20. * units.min
		tau = 21. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		# Four generations means one event from first generation, two from the
		# second, 4 from the third, and 8 from the fourth - total of 15 events, 2 forks per event = 30 replisomes
		self.assertTrue(len(replisome_state["coordinates"]) == len(replisome_state["right_replichore"]) == len(replisome_state["domain_index"]) == 30)


	@noseAttrib.attr('replicationTest')
	def test_fork_coordinates(self):

		C = 70. * units.min
		D = 25. * units.min
		tau = 30. * units.min
		replichore_length = 2319838 * units.nt

		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)
		fork_coordinates = replisome_state["coordinates"]

		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.convertNoUnitToNumber(ratio)
			self.assertTrue(fork_coordinates[2**n - 2] == -fork_coordinates[2**n - 1] == np.floor(ratio*(replichore_length.asNumber())))
			n += 1

		# Length of DNA to be replicated is 1
		C = 70. * units.min
		D = 25. * units.min
		tau = 30. * units.min
		replichore_length = 1 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)
		fork_coordinates = replisome_state["coordinates"]
	
		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.convertNoUnitToNumber(ratio)
			self.assertTrue(fork_coordinates[2**n - 2] == -fork_coordinates[2**n - 1] == np.floor(ratio*(replichore_length.asNumber())))
			n += 1

		# Length of DNA to be replicated is 0
		C = 70. * units.min
		D = 25. * units.min
		tau = 30. * units.min
		replichore_length = 0 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)
		fork_coordinates = replisome_state["coordinates"]

		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.convertNoUnitToNumber(ratio)
			self.assertTrue(fork_coordinates[2**n - 2] == -fork_coordinates[2**n - 1] == np.floor(ratio*(replichore_length.asNumber())))
			n += 1

		# Length of DNA to be replicated is not a whole number
		C = 70. * units.min
		D = 25. * units.min
		tau = 30. * units.min
		replichore_length = 2319838.3 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)
		fork_coordinates = replisome_state["coordinates"]

		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.convertNoUnitToNumber(ratio)
			self.assertTrue(fork_coordinates[2**n - 2] == -fork_coordinates[2**n - 1] == np.floor(ratio*(replichore_length.asNumber())))
			n += 1


	@noseAttrib.attr('replicationTest')
	def test_maximum_domain_index(self):

		# When (C + D) / tau is less than one, no replication will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 90. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		self.assertEqual(domain_state["domain_index"].max(), 0)


		# When (C + D) / tau is between one and two, one replication generation will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 60. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		self.assertEqual(domain_state["domain_index"].max(), 2)


		# When (C + D) / tau is between two and three, two replication generations will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 30. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		# There should be equal numbers of zeros and ones in the chromosomeIndex array (excepting the first four)
		self.assertEqual(domain_state["domain_index"].max(), 6)

		
		# When (C + D) / tau is between three and four, three replication generations will have started
		C = 50. * units.min
		D = 19. * units.min
		tau = 20. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		# There should be equal numbers of zeros and ones in the chromosomeIndex array (excepting the first four)
		self.assertEqual(domain_state["domain_index"].max(), 14)


		# When (C + D) / tau is between four and five, four replication generations will have started
		C = 70. * units.min
		D = 20. * units.min
		tau = 21. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		# There should be equal numbers of zeros and ones in the chromosomeIndex array (excepting the first four)
		self.assertEqual(domain_state["domain_index"].max(), 30)


	@noseAttrib.attr('replicationTest')
	def test_determine_chromosome_state_inputs(self):

		# The D period must be shorter than tau
		with self.assertRaises(AssertionError) as context:
			C = 50. * units.min
			D = 20. * units.min
			tau = 19. * units.min
			replichore_length = 2319838  * units.nt
			oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)
		self.assertEqual(context.exception.message, "The D period must be shorter than the doubling time tau.")


		# No inputs can have negative values
		with self.assertRaises(AssertionError) as context:
			C = -50. * units.min
			D = 20. * units.min
			tau = 60. * units.min
			replichore_length = 2319838  * units.nt
			oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)
		self.assertEqual(context.exception.message, "C value can't be negative.")

		with self.assertRaises(AssertionError) as context:
			C = 50. * units.min
			D = -20. * units.min
			tau = 60. * units.min
			replichore_length = 2319838  * units.nt
			oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)
		self.assertEqual(context.exception.message, "D value can't be negative.")

		with self.assertRaises(AssertionError) as context:
			C = 50. * units.min
			D = 20. * units.min
			tau = -60. * units.min
			replichore_length = 2319838  * units.nt
			oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)
		self.assertEqual(context.exception.message, "tau value can't be negative.")

		with self.assertRaises(AssertionError) as context:
			C = 50. * units.min
			D = 20. * units.min
			tau = 60. * units.min
			replichore_length = -2319838  * units.nt
			oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)
		self.assertEqual(context.exception.message, "replichore_length value can't be negative.")


	@noseAttrib.attr('replicationTest')
	def test_num_oriCs(self):

		# When (C + D) / tau is less than one, no replication will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 100. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		self.assertEqual(len(oric_state["domain_index"]), 1)

		# When (C + D) / tau is less than one, no replication will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 90. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		self.assertEqual(len(oric_state["domain_index"]), 1)

		# When (C + D) / tau is between one and two, one replication generation will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 60. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		self.assertEqual(len(oric_state["domain_index"]), 2)

		# When (C + D) / tau is between two and three, two replication generations will have started
		C = 50. * units.min
		D = 27. * units.min
		tau = 30. * units.min
		replichore_length = 2319838 * units.nt
		oric_state, replisome_state, domain_state = determine_chromosome_state(C, D, tau, replichore_length, -1)

		# Two generations means one event from first generation and two events
		# from the second - total of 3 events, plus one oriC for the original
		# chromosome = 4 total oriC
		self.assertEqual(len(oric_state["domain_index"]), 4)


if __name__ == '__main__':
	unittest.main()
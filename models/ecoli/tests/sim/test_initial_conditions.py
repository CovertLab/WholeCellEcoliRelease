"""
Test_InitialConditions.py
Tests the initial conditions code of the model.

@author: Morgan Paull
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

from  models.ecoli.sim.initial_conditions import determineChromosomeState, determineNumOriC


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
	def test_num_DNA_polys(self):

		# When (C + D) / tau is less than one, no replication will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 90 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

		assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(chromosomeIndex) == 0)

		# When (C + D) / tau is between one and two, one replication generation will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 60 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

		assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(chromosomeIndex) == 4)

		# When (C + D) / tau is between two and three, two replication generations will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 30 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

		# Two generations means one event from first generation and two events
		# from the second - total of 3 events, 4 DNA polys per event = 12 DNA polys
		assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(chromosomeIndex) == 12)

		# When (C + D) / tau is between three and four, three replication generations will have started
		C = 50 * units.min
		D = 19 * units.min
		tau = 20 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

		# Three generations means one event from first generation, two from the
		# second, and 4 from the third - total of 7 events, 4 DNA polys per event = 28 DNA polys
		assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(chromosomeIndex) == 28)

		# When (C + D) / tau is between four and five, four replication generations will have started
		C = 70 * units.min
		D = 20 * units.min
		tau = 21 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

		# Four generations means one event from first generation, two from the
		# second, 4 from the third, and 8 from the fourth - total of 15 events, 4 DNA polys per event = 60 DNA polys
		assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(chromosomeIndex) == 60)





	@noseAttrib.attr('replicationTest')
	def test_DNA_polys_replicated_lengths(self):

		C = 70 * units.min
		D = 25 * units.min
		tau = 30 * units.min
		replication_length = 2319838 * units.nt

		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1;
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.convertNoUnitToNumber(ratio)
			assert ( sequenceLength[4*(n-1)] == sequenceLength[4*(n-1)+1] == sequenceLength[4*(n-1)+2] == sequenceLength[4*(n-1)+3] == np.floor(ratio*(replication_length.asNumber())) )
			n += 1

		# Length of DNA to be replicated is 1
		C = 70 * units.min
		D = 25 * units.min
		tau = 30 * units.min
		replication_length = 1 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
	
		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1;
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.convertNoUnitToNumber(ratio)
			assert ( sequenceLength[4*(n-1)] == sequenceLength[4*(n-1)+1] == sequenceLength[4*(n-1)+2] == sequenceLength[4*(n-1)+3] == np.floor(ratio*(replication_length.asNumber())) )
			n += 1

		# Length of DNA to be replicated is 0
		C = 70 * units.min
		D = 25 * units.min
		tau = 30 * units.min
		replication_length = 0 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
	
		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1;
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.convertNoUnitToNumber(ratio)
			assert ( sequenceLength[4*(n-1)] == sequenceLength[4*(n-1)+1] == sequenceLength[4*(n-1)+2] == sequenceLength[4*(n-1)+3] == np.floor(ratio*(replication_length.asNumber())) )
			n += 1

		# Length of DNA to be replicated is infinite
		C = 70 * units.min
		D = 25 * units.min
		tau = 30 * units.min
		replication_length = np.inf * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1;
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.convertNoUnitToNumber(ratio)
			assert ( sequenceLength[4*(n-1)] == sequenceLength[4*(n-1)+1] == sequenceLength[4*(n-1)+2] == sequenceLength[4*(n-1)+3] == np.floor(ratio*(replication_length.asNumber())) )
			n += 1

		# Length of DNA to be replicated is not a whole number
		C = 70 * units.min
		D = 25 * units.min
		tau = 30 * units.min
		replication_length = 2319838.3 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

		limit = np.floor((C.asNumber() + D.asNumber())/tau.asNumber())
		n = 1;
		while n <= limit:
			ratio = (1 - ((n*tau - D)/(C)))
			ratio = units.convertNoUnitToNumber(ratio)
			assert ( sequenceLength[4*(n-1)] == sequenceLength[4*(n-1)+1] == sequenceLength[4*(n-1)+2] == sequenceLength[4*(n-1)+3] == np.floor(ratio*(replication_length.asNumber())) )
			n += 1



	@noseAttrib.attr('replicationTest')
	def test_DNA_polys_division_characteristics(self):

		# When (C + D) / tau is less than one, no replication will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 90 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
		
		# Array should be of length zero if no replication events have started
		assert( sum(chromosomeIndex) == len(chromosomeIndex) == 0)


		# When (C + D) / tau is between one and two, one replication generation will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 60 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

		assert( sum(chromosomeIndex) == (len(chromosomeIndex)-4)/2 )


		# When (C + D) / tau is between two and three, two replication generations will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 30 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

		# There should be equal numbers of zeros and ones in the chromosomeIndex array (excepting the first four)
		assert( sum(chromosomeIndex) == (len(chromosomeIndex)-4)/2 )

		
		# When (C + D) / tau is between three and four, three replication generations will have started
		C = 50 * units.min
		D = 19 * units.min
		tau = 20 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

		# There should be equal numbers of zeros and ones in the chromosomeIndex array (excepting the first four)
		assert( sum(chromosomeIndex) == (len(chromosomeIndex)-4)/2 )


		# When (C + D) / tau is between four and five, four replication generations will have started
		C = 70 * units.min
		D = 20 * units.min
		tau = 21 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)

		# There should be equal numbers of zeros and ones in the chromosomeIndex array (excepting the first four)
		assert( sum(chromosomeIndex) == (len(chromosomeIndex)-4)/2 )


	@noseAttrib.attr('replicationTest')
	def test_determineChromosomeState_inputs(self):

		# The D period must be shorter than tau
		with self.assertRaises(AssertionError) as context:
			C = 50 * units.min
			D = 20  * units.min
			tau = 19 * units.min
			replication_length = 2319838  * units.nt
			sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
		self.assertEqual(context.exception.message, "The D period must be shorter than the doubling time tau.")


		# No inputs can have negative values
		with self.assertRaises(AssertionError) as context:
			C = -50 * units.min
			D = 20  * units.min
			tau = 60 * units.min
			replication_length = 2319838  * units.nt
			sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
		self.assertEqual(context.exception.message, "C value can't be negative.")

		with self.assertRaises(AssertionError) as context:
			C = 50 * units.min
			D = -20  * units.min
			tau = 60 * units.min
			replication_length = 2319838  * units.nt
			sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
		self.assertEqual(context.exception.message, "D value can't be negative.")

		with self.assertRaises(AssertionError) as context:
			C = 50 * units.min
			D = 20  * units.min
			tau = -60 * units.min
			replication_length = 2319838  * units.nt
			sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
		self.assertEqual(context.exception.message, "tau value can't be negative.")

		with self.assertRaises(AssertionError) as context:
			C = 50 * units.min
			D = 20  * units.min
			tau = 60 * units.min
			replication_length = -2319838  * units.nt
			sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
		self.assertEqual(context.exception.message, "replication_length value can't be negative.")


	@noseAttrib.attr('replicationTest')
	def test_determineChromosomeState_input_units(self):

		# Inputs with units (should not throw an error)
		C = 50 * units.min
		D = 27 * units.min
		tau = 90 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)


		# Inputs with no units (expect an error)
		with self.assertRaises(AssertionError) as context:
			C = 50
			D = 27  * units.min
			tau = 90  * units.min
			replication_length = 2319838  * units.nt
			sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
		self.assertEqual(context.exception.message, "C must have units")

		with self.assertRaises(AssertionError) as context:
			C = 50 * units.min
			D = 27
			tau = 90  * units.min
			replication_length = 2319838  * units.nt
			sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
		self.assertEqual(context.exception.message, "D must have units")

		with self.assertRaises(AssertionError) as context:
			C = 50 * units.min
			D = 27  * units.min
			tau = 90
			replication_length = 2319838  * units.nt
			sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
		self.assertEqual(context.exception.message, "tau must have units")

		with self.assertRaises(AssertionError) as context:
			C = 50 * units.min
			D = 27  * units.min
			tau = 90  * units.min
			replication_length = 2319838
			sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
		self.assertEqual(context.exception.message, "replication_length must have units of units.nt.")


		# Inputs with incorrect units (expect an error)
		with self.assertRaises(AssertionError) as context:
			C = 50 * units.min
			D = 27  * units.min
			tau = 90  * units.min
			replication_length = 2319838  * units.min
			sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
		self.assertEqual(context.exception.message, "replication_length must have units of units.nt.")

		with self.assertRaises(AssertionError) as context:
			C = 50 * units.min
			D = 27  * units.min
			tau = 90  * units.min
			replication_length = 2319838  * units.s
			sequenceIdx, sequenceLength, replicationRound, chromosomeIndex = determineChromosomeState(C, D, tau, replication_length)
		self.assertEqual(context.exception.message, "replication_length must have units of units.nt.")


	@noseAttrib.attr('replicationTest')
	def test_num_oriCs(self):

		# When (C + D) / tau is less than one, no replication will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 100 * units.min
		numOric = determineNumOriC(C, D, tau)

		assert(numOric == 1)

		# When (C + D) / tau is less than one, no replication will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 90 * units.min
		numOric = determineNumOriC(C, D, tau)

		assert(numOric == 1)

		# When (C + D) / tau is between one and two, one replication generation will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 60 * units.min
		numOric = determineNumOriC(C, D, tau)

		assert(numOric == 2)

		# When (C + D) / tau is between two and three, two replication generations will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 30 * units.min
		numOric = determineNumOriC(C, D, tau)

		# Two generations means one event from first generation and two events
		# from the second - total of 3 events, plus one oriC for the original
		# chromosome = 4 total oriC
		assert(numOric == 4)

		# When (C + D) / tau is between three and four, three replication generations will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 20 * units.min
		numOric = determineNumOriC(C, D, tau)

		# Three generations means one event from first generation, two from the
		# second, and 4 from the third - total of 7 events, plus one oriC for 
		# the original chromosome = 8 total oriC
		assert(numOric == 8)

		# When (C + D) / tau is between four and five, four replication generations will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 17 * units.min
		numOric = determineNumOriC(C, D, tau)

		# Four generations means one event from first generation, two from the
		# second, 4 from the third, and 8 from the fourth - total of 15 events,
		# plus one oriC for the original chromosome = 16 total oriC
		assert(numOric == 16)

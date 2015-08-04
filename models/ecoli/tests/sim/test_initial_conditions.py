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

from  models.ecoli.sim.initial_conditions.py import determineChromosomeState


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
		sequenceIdx, sequenceLength, replicationRound, replicationDivision, numOric = determineChromosomeState(C, D, tau, replication_length)

		assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(replicationDivision) == 0)

		# When (C + D) / tau is between one and two, one replication generation will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 60 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, replicationDivision, numOric = determineChromosomeState(C, D, tau, replication_length)

		assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(replicationDivision) == 4)

		# When (C + D) / tau is between two and three, two replication generations will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 30 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, replicationDivision, numOric = determineChromosomeState(C, D, tau, replication_length)

		# Two generations means one event from first generation and two events
		# from the second - total of 3 events, 4 DNA polys per event = 12 DNA polys
		assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(replicationDivision) == 12)

		# When (C + D) / tau is between three and four, three replication generations will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 20 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, replicationDivision, numOric = determineChromosomeState(C, D, tau, replication_length)

		# Three generations means one event from first generation, two from the
		# second, and 4 from the third - total of 7 events, 4 DNA polys per event = 28 DNA polys
		assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(replicationDivision) == 28)

		# When (C + D) / tau is between four and five, four replication generations will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 17 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, replicationDivision, numOric = determineChromosomeState(C, D, tau, replication_length)

		# Four generations means one event from first generation, two from the
		# second, 4 from the third, and 8 from the fourth - total of 15 events, 4 DNA polys per event = 60 DNA polys
		assert(len(sequenceIdx) == len(sequenceLength) == len(replicationRound) == len(replicationDivision) == 60)


	@noseAttrib.attr('replicationTest')
	def test_num_oriCs(self):

		# When (C + D) / tau is less than one, no replication will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 90 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, replicationDivision, numOric = determineChromosomeState(C, D, tau, replication_length)

		assert(len(numOric) == 1)

		# When (C + D) / tau is between one and two, one replication generation will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 60 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, replicationDivision, numOric = determineChromosomeState(C, D, tau, replication_length)

		assert(len(numOric) == 2)

		# When (C + D) / tau is between two and three, two replication generations will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 30 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, replicationDivision, numOric = determineChromosomeState(C, D, tau, replication_length)

		# Two generations means one event from first generation and two events
		# from the second - total of 3 events, plus one oriC for the original
		# chromosome = 4 total oriC
		assert(len(numOric) == 4)

		# When (C + D) / tau is between three and four, three replication generations will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 20 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, replicationDivision, numOric = determineChromosomeState(C, D, tau, replication_length)

		# Three generations means one event from first generation, two from the
		# second, and 4 from the third - total of 7 events, plus one oriC for 
		# the original chromosome = 8 total oriC
		assert(len(numOric) == 8)

		# When (C + D) / tau is between four and five, four replication generations will have started
		C = 50 * units.min
		D = 27 * units.min
		tau = 17 * units.min
		replication_length = 2319838 * units.nt
		sequenceIdx, sequenceLength, replicationRound, replicationDivision, numOric = determineChromosomeState(C, D, tau, replication_length)

		# Four generations means one event from first generation, two from the
		# second, 4 from the third, and 8 from the fourth - total of 15 events,
		# plus one oriC for the original chromosome = 16 total oriC
		assert(len(numOric) == 16)
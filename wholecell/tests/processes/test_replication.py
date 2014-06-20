#!/usr/bin/env python

"""
Test replication.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Chemical Engineering, Stanford University
@date: Created 4/6/2014
"""
import wholecell.processes.replication

import numpy as np

import nose.plugins.attrib as noseAttrib
import nose.tools as noseTools
import unittest

class Test_replication(unittest.TestCase):
	
	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def setUp(self):
										# Initial set up
										# - = fork initial position
										# x = terC location

		                          		######-#-#########x####
		self.genomeSequence = np.array([0,1,2,3,1,3,2,1,2,3,1,2], dtype=np.int8)
		self.dnaPolymeraseElongationRate = 2
		self.genomeLength = len(self.genomeSequence)
		self.nPolymerase = 4
		self.allChromosomeLocation = np.array([3,3,4,4])
		self.allDirectionIsPositive = np.array([False, False, True, True])
		self.allIsLeading = np.array([True, False, True, False])
		self.tercCenter = self.genomeLength - 2

	def tearDown(self):
		pass

	@noseAttrib.attr('replicationTest')
	@noseAttrib.attr('smalltest')
	def test_calculateUpcomingSequence(self):
		# Test normal moves not encountering the terC
		chromosomeLocation = [0,0,5,5]
		directionIsPositive = [True,True,False,False]
		isLeading = [True,False,True,False]
		elongationRate = [5,5,6,6]

		comparison_sequence = [
			self.genomeSequence[0:5],
			wholecell.processes.replication.reverseComplement(self.genomeSequence[0:5]),
			self.genomeSequence[np.arange(5,-1,-1) % self.genomeLength],
			wholecell.processes.replication.reverseComplement(self.genomeSequence[np.arange(5,-1,-1) % self.genomeLength]),
			]

		for i in range(len(chromosomeLocation)):
			upcomingSequence = wholecell.processes.replication.calculateUpcomingSequence(
				chromosomeLocation = chromosomeLocation[i],
				directionIsPositive = directionIsPositive[i],
				isLeading = isLeading[i],
				elongationRate = elongationRate[i],
				genomeLength = self.genomeLength,
				genomeSequence = self.genomeSequence,
				tercCenter = self.tercCenter
				)

			self.assertEqual(upcomingSequence.tolist(), comparison_sequence[i].tolist())

		# Test encounter with terC
		chromosomeLocation = [0,self.genomeLength - 1]
		directionIsPositive = [True,False]
		isLeading = [True,True]
		elongationRate = [self.genomeLength,self.genomeLength]

		comparison_sequence = [
			self.genomeSequence[0:self.tercCenter],
			self.genomeSequence[self.genomeLength - 1:self.tercCenter:-1],
			]

		for i in range(len(chromosomeLocation)):
			upcomingSequence = wholecell.processes.replication.calculateUpcomingSequence(
				chromosomeLocation = chromosomeLocation[i],
				directionIsPositive = directionIsPositive[i],
				isLeading = isLeading[i],
				elongationRate = elongationRate[i],
				genomeLength = self.genomeLength,
				genomeSequence = self.genomeSequence,
				tercCenter = self.tercCenter
				)

			self.assertEqual(upcomingSequence.tolist(), comparison_sequence[i].tolist())

	@noseAttrib.attr('replicationTest')
	@noseAttrib.attr('smalltest')
	def	test_reverseComplement(self):
		self.assertEqual(
			wholecell.processes.replication.reverseComplement(self.genomeSequence).tolist(),
			[3,2,1,0,2,0,1,2,1,0,2,1]
			)

	@noseAttrib.attr('replicationTest')
	@noseAttrib.attr('smalltest')
	def	test_buildSequenceMatrix(self):
		return
		test_sequenceMatrix = wholecell.processes.replication.buildSequenceMatrix(
			self.nPolymerase,
			self.allChromosomeLocation,
			self.allDirectionIsPositive,
			self.allIsLeading,
			self.dnaPolymeraseElongationRate,
			self.genomeLength,
			self.genomeSequence,
			self.tercCenter
			)

		sequenceMatrix = np.empty((self.nPolymerase, self.dnaPolymeraseElongationRate), np.int8)
		sequenceMatrix.fill(wholecell.processes.replication.PAD_VALUE)
		sequenceMatrix[0,:] = [2,1]
		sequenceMatrix[1,:] = [1,2]
		sequenceMatrix[2,:] = [3,2]
		sequenceMatrix[3,:] = [0,1]
		import ipdb; ipdb.set_trace()



	@noseAttrib.attr('replicationTest')
	@noseAttrib.attr('smalltest')
	def test_calculatePolymerasePositionUpdate(self):
		# Positive direction
		currentPosition = 0
		directionIsPositive = True
		difference = 4

		newPos = wholecell.processes.replication.calculatePolymerasePositionUpdate(
			currentPosition, directionIsPositive, difference, self.genomeLength
			)

		self.assertEqual(newPos, 4)

		# Negative direction
		currentPosition = len(self.genomeSequence) - 1
		directionIsPositive = False
		difference = 4

		newPos = wholecell.processes.replication.calculatePolymerasePositionUpdate(
			currentPosition, directionIsPositive, difference, self.genomeLength
			)

		self.assertEqual(newPos, len(self.genomeSequence) - 1 - 4)

		# Positive loop condition
		currentPosition = len(self.genomeSequence) - 2
		directionIsPositive = True
		difference = 4

		newPos = wholecell.processes.replication.calculatePolymerasePositionUpdate(
			currentPosition, directionIsPositive, difference, self.genomeLength
			)

		self.assertEqual(newPos, 2)

		# Negative loop condition
		currentPosition = 2
		directionIsPositive = False
		difference = 4

		newPos = wholecell.processes.replication.calculatePolymerasePositionUpdate(
			currentPosition, directionIsPositive, difference, self.genomeLength
			)

		self.assertEqual(newPos, len(self.genomeSequence) - 1 - 1)
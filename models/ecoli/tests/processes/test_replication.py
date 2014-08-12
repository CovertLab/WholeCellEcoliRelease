#!/usr/bin/env python

"""
Test replication.py

@author: Nick Ruggero
@organization: Covert Lab, Department of Chemical Engineering, Stanford University
@date: Created 4/6/2014
"""
import models.ecoli.processes.replication

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

		geneEndCoordinate = np.array([3,10])
		self.bufferedGeneEndCoordinate = np.concatenate(
			[geneEndCoordinate - self.genomeLength, geneEndCoordinate, geneEndCoordinate + self.genomeLength]
			)

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
			models.ecoli.processes.replication.reverseComplement(self.genomeSequence[0:5]),
			self.genomeSequence[np.arange(5,-1,-1) % self.genomeLength],
			models.ecoli.processes.replication.reverseComplement(self.genomeSequence[np.arange(5,-1,-1) % self.genomeLength]),
			]

		for i in range(len(chromosomeLocation)):
			upcomingSequence = models.ecoli.processes.replication.calculateUpcomingSequence(
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
			upcomingSequence = models.ecoli.processes.replication.calculateUpcomingSequence(
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
			models.ecoli.processes.replication.reverseComplement(self.genomeSequence).tolist(),
			[3,2,1,0,2,0,1,2,1,0,2,1]
			)

	@noseAttrib.attr('replicationTest')
	@noseAttrib.attr('smalltest')
	def	test_buildSequenceMatrix(self):
		test_sequenceMatrix = models.ecoli.processes.replication.buildSequenceMatrix(
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
		sequenceMatrix.fill(models.ecoli.processes.replication.PAD_VALUE)
		sequenceMatrix[0,:] = [3,2]
		sequenceMatrix[1,:] = [0,1]
		sequenceMatrix[2,:] = [1,3]
		sequenceMatrix[3,:] = [2,0]

		self.assertTrue(np.all(test_sequenceMatrix == sequenceMatrix))

	@noseAttrib.attr('replicationTest')
	@noseAttrib.attr('smalltest')
	def test_calculatePolymerasePositionUpdate(self):
		# Positive direction
		currentPosition = 0
		directionIsPositive = True
		difference = 4

		newPos = models.ecoli.processes.replication.calculatePolymerasePositionUpdate(
			currentPosition, directionIsPositive, difference, self.genomeLength
			)

		self.assertEqual(newPos, 4)

		# Negative direction
		currentPosition = len(self.genomeSequence) - 1
		directionIsPositive = False
		difference = 4

		newPos = models.ecoli.processes.replication.calculatePolymerasePositionUpdate(
			currentPosition, directionIsPositive, difference, self.genomeLength
			)

		self.assertEqual(newPos, len(self.genomeSequence) - 1 - 4)

		# Positive loop condition
		currentPosition = len(self.genomeSequence) - 2
		directionIsPositive = True
		difference = 4

		newPos = models.ecoli.processes.replication.calculatePolymerasePositionUpdate(
			currentPosition, directionIsPositive, difference, self.genomeLength
			)

		self.assertEqual(newPos, 2)

		# Negative loop condition
		currentPosition = 2
		directionIsPositive = False
		difference = 4

		newPos = models.ecoli.processes.replication.calculatePolymerasePositionUpdate(
			currentPosition, directionIsPositive, difference, self.genomeLength
			)

		self.assertEqual(newPos, len(self.genomeSequence) - 1 - 1)

	@noseAttrib.attr('replicationTest')
	@noseAttrib.attr('smalltest')
	def test_calculateReplicatedGenes(self):
		replicatedGenes = models.ecoli.processes.replication.calculateReplicatedGenes(
			currentPosition = 0,
			directionIsPositive = True,
			difference = 5,
			bufferedGeneEndCoordinate = self.bufferedGeneEndCoordinate
			)
		self.assertTrue(np.all(replicatedGenes == np.array([True, False])))

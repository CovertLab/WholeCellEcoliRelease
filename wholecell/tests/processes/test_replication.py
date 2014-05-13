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
		self.sequence = np.array(list('hello_world_i_am_genome'))
		self.elongationRate = 5
		self.genomeLength = len(self.sequence)

	def tearDown(self):
		pass

	@noseAttrib.attr('replicationTest')
	@noseAttrib.attr('smalltest')
	def test_calculateSequence(self):
		# Positive direction
		chromosomeLocation = 0
		directionIsPositive = True
		
		seq = wholecell.processes.replication.calculateSequence(
			chromosomeLocation, directionIsPositive, self.elongationRate, self.sequence, self.genomeLength
			)

		self.assertEqual(seq.tolist(), list('hello'))
		
		# Negative direction
		chromosomeLocation = len(self.sequence) - 1
		directionIsPositive = False
		
		seq = wholecell.processes.replication.calculateSequence(
			chromosomeLocation, directionIsPositive, self.elongationRate, self.sequence, self.genomeLength
			)
		self.assertEqual(seq.tolist(), list('enome'[::-1]))
		
		# Loop condition positive direction
		chromosomeLocation = len(self.sequence) - 2
		directionIsPositive = True

		seq = wholecell.processes.replication.calculateSequence(
			chromosomeLocation, directionIsPositive, self.elongationRate, self.sequence, self.genomeLength
			)
		self.assertEqual(seq.tolist(), list('me' + 'hel'))
		
		# Loop condition negative direction
		chromosomeLocation = 2
		directionIsPositive = False

		seq = wholecell.processes.replication.calculateSequence(
			chromosomeLocation, directionIsPositive, self.elongationRate, self.sequence, self.genomeLength
			)

		self.assertEqual(seq.tolist(), list('hel'[::-1] + 'me'[::-1]))

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
		currentPosition = len(self.sequence) - 1
		directionIsPositive = False
		difference = 4

		newPos = wholecell.processes.replication.calculatePolymerasePositionUpdate(
			currentPosition, directionIsPositive, difference, self.genomeLength
			)

		self.assertEqual(newPos, len(self.sequence) - 1 - 4)

		# Positive loop condition
		currentPosition = len(self.sequence) - 2
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

		self.assertEqual(newPos, len(self.sequence) - 1 - 1)
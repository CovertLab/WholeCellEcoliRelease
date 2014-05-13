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
		self.sequence = 'hello_world_i_am_genome'
		self.elongationRate = 5
		self.genomeLength = len(self.sequence)

	def tearDown(self):
		pass

	@noseAttrib.attr('replicationTest')
	@noseAttrib.attr('smalltest')
	def test_calculateSequence(self):
		# Positive direction
		chromosomeLocation = 1
		directionIsPositive = True
		
		seq = wholecell.processes.replication.calculateSequence(
			chromosomeLocation, directionIsPositive, self.elongationRate, self.sequence, self.genomeLength
			)
		self.assertEqual(seq, 'hello')

		# Negative direction
		chromosomeLocation = len(self.sequence)
		directionIsPositive = False
		
		seq = wholecell.processes.replication.calculateSequence(
			chromosomeLocation, directionIsPositive, self.elongationRate, self.sequence, self.genomeLength
			)
		self.assertEqual(seq, 'enome'[::-1])
		
		# Loop condition positive direction
		chromosomeLocation = len(self.sequence) - 2
		directionIsPositive = True

		seq = wholecell.processes.replication.calculateSequence(
			chromosomeLocation, directionIsPositive, self.elongationRate, self.sequence, self.genomeLength
			)
		self.assertEqual(seq, 'me' + 'hel')
		return
		# Loop condition negative direction
		chromosomeLocation = 2
		directionIsPositive = False

		seq = wholecell.processes.replication.calculateSequence(
			chromosomeLocation, directionIsPositive, self.elongationRate, self.sequence, self.genomeLength
			)

		self.assertEqual(seq, 'AAC'[::-1] + 'CA'[::-1])
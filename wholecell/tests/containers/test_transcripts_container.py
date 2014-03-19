'''
test_chromosome_container.py

Tests for the ChromosomeContainer class.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 3/17/2014
'''

from __future__ import division

import unittest

import numpy as np
import nose.plugins.attrib as noseAttrib

from wholecell.containers.transcripts_container import TranscriptsContainer

ARRAY_LENGTH = 1000

MOLECULE_ATTRIBUTES = {}

class Test_TranscriptsContainer(unittest.TestCase):
	@classmethod
	def setupClass(cls):
		pass


	@classmethod
	def tearDownClass(cls):
		pass


	def setUp(self):
		self.container = createContainer()


	def tearDown(self):
		pass


	@noseAttrib.attr('smalltest', 'transcripts', 'containerObject')
	def test_transcript_new(self):
		size = 100

		transcript = self.container.transcriptNew(0, size)

		transcriptPosition = np.where(
			self.container._array == (transcript.attr('_globalIndex') + self.container._offset)
			)[0]

		self.assertEqual(transcriptPosition.size, 1)

		self.assertEqual(
			transcript.attr('_transcriptPosition'),
			transcriptPosition[0]
			)

		reservedPositions = np.where(
			self.container._array == self.container._reserved
			)[0]

		self.assertEqual(reservedPositions.size, size)


	@noseAttrib.attr('smalltest', 'transcripts', 'containerObject')
	def test_transcript_extend(self):
		transcript = self.container.transcriptNew(0, 100)

		self.container.transcriptExtend(transcript, 20)

		positions = np.where(self.container._array == self.container._empty)[0]

		self.assertEqual(positions.size, 20)

		self.container.transcriptExtend(transcript, 20)

		positions = np.where(self.container._array == self.container._empty)[0]

		self.assertEqual(positions.size, 40)


	@noseAttrib.attr('smalltest', 'transcripts', 'containerObject')
	def test_transcript_del(self):
		transcript = self.container.transcriptNew(0, 100)

		self.container.transcriptExtend(transcript, 20)

		self.container.transcriptDel(transcript)

		positions = np.where(self.container._array == self.container._empty)[0]

		self.assertEqual(positions.size, 0)

		self.assertEqual(
			self.container,
			createContainer()
			)


	@noseAttrib.attr('smalltest', 'transcripts', 'containerObject')
	def test_many_transcripts(self):
		nTranscripts = 25
		size = 10

		transcripts = [self.container.transcriptNew(0, size)
			for i in xrange(nTranscripts)]

		for transcript in transcripts:
			transcriptPosition = np.where(
				self.container._array == (transcript.attr('_globalIndex') + self.container._offset)
				)[0]

			self.assertEqual(transcriptPosition.size, 1)

			self.assertEqual(
				transcript.attr('_transcriptPosition'),
				transcriptPosition[0]
				)

		reservedPositions = np.where(
			self.container._array == self.container._reserved
			)[0]

		self.assertEqual(reservedPositions.size, nTranscripts * size)
		
		extend = 20

		for transcript in transcripts:
			self.container.transcriptExtend(transcript, extend)

		openPositions = np.where(
			self.container._array == self.container._empty
			)[0]

		self.assertEqual(openPositions.size, nTranscripts * extend)

		for transcript in transcripts:
			self.container.transcriptDel(transcript)

		self.assertEqual(
			self.container,
			createContainer()
			)


def createContainer():
	container = TranscriptsContainer(ARRAY_LENGTH, MOLECULE_ATTRIBUTES)

	return container

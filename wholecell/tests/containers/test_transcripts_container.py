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

MOLECULE_ATTRIBUTES = {
	'Ribosome':{},
	'RNAse':{}
	}

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
			transcript.attr('_transPosition'),
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
				transcript.attr('_transPosition'),
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


	@noseAttrib.attr('smalltest', 'transcripts', 'containerObject')
	def test_bind_molecule(self):
		transcript = self.container.transcriptNew(0, 100)
		self.container.transcriptExtend(transcript, 100)

		molecule = self.container.moleculeNew('Ribosome')

		self.container.moleculeLocationIs(molecule, transcript, 50, '+',
			10, 5)

		moleculePosition = np.where(
			self.container._array == (molecule.attr('_globalIndex') + self.container._offset)
			)[0]

		relativePosition = moleculePosition - (transcript.attr('_transPosition') + 1)

		self.assertEqual(
			relativePosition.tolist(),
			[45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59]
			)


	@noseAttrib.attr('smalltest', 'transcripts', 'containerObject')
	def test_bind_molecule_reverse(self):
		transcript = self.container.transcriptNew(0, 100)
		self.container.transcriptExtend(transcript, 100)

		molecule = self.container.moleculeNew('Ribosome')

		self.container.moleculeLocationIs(molecule, transcript, 50, '-',
			10, 5)

		moleculePosition = np.where(
			self.container._array == (molecule.attr('_globalIndex') + self.container._offset)
			)[0]

		relativePosition = moleculePosition - (transcript.attr('_transPosition') + 1)

		self.assertEqual(
			relativePosition.tolist(),
			[41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55]
			)


	@noseAttrib.attr('smalltest', 'transcripts', 'containerObject')
	def test_move_molecule(self):
		transcript = self.container.transcriptNew(0, 100)
		self.container.transcriptExtend(transcript, 100)

		molecule = self.container.moleculeNew('Ribosome')

		self.container.moleculeLocationIs(molecule, transcript, 45, '+',
			10, 5)

		self.container.moleculeLocationIs(molecule, transcript, 50, '+',
			10, 5)

		moleculePosition = np.where(
			self.container._array == (molecule.attr('_globalIndex') + self.container._offset)
			)[0]

		relativePosition = moleculePosition - (transcript.attr('_transPosition') + 1)

		self.assertEqual(
			relativePosition.tolist(),
			[45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59]
			)


	@noseAttrib.attr('smalltest', 'transcripts', 'containerObject')
	def test_many_transcripts_many_molecules(self):
		nPairs = 25
		sizeAllocated = 0
		extend = 20

		pairs = []

		for i in xrange(nPairs):
			transcript = self.container.transcriptNew(0, sizeAllocated)
			self.container.transcriptExtend(transcript, extend)

			molecule = self.container.moleculeNew('Ribosome')

			self.container.moleculeLocationIs(molecule, transcript, 10, '+',
				10, 5)

			pairs.append((transcript, molecule))

		for transcript, molecule in pairs:
			moleculePosition = np.where(
				self.container._array == (molecule.attr('_globalIndex') + self.container._offset)
				)[0]

			relativePosition = moleculePosition - (transcript.attr('_transPosition') + 1)
			
			self.assertEqual(
				relativePosition.tolist(),
				[5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
				)


	@noseAttrib.attr('smalltest', 'transcripts', 'containerObject')
	def test_moleculesBound(self):
		transcript = self.container.transcriptNew(0, 100)
		self.container.transcriptExtend(transcript, 100)

		positions = [10, 30, 50]
		molecules = []

		for position in positions:
			molecule = self.container.moleculeNew('Ribosome')

			self.container.moleculeLocationIs(molecule, transcript, position, '+',
				10, 5)

			molecules.append(molecule)

		unboundMolecule = self.container.moleculeNew('Ribosome')

		self.assertEqual(
			set(self.container.moleculesBound()),
			set(molecules)
			)


	@noseAttrib.attr('smalltest', 'transcripts', 'containerObject')
	def test_moleculesUnbound(self):
		transcript = self.container.transcriptNew(0, 100)
		self.container.transcriptExtend(transcript, 100)

		positions = [10, 30, 50]
		molecules = []

		for position in positions:
			molecule = self.container.moleculeNew('Ribosome')

			self.container.moleculeLocationIs(molecule, transcript, position, '+',
				10, 5)

			molecules.append(molecule)

		unboundMolecule = self.container.moleculeNew('Ribosome')

		self.assertEqual(
			set(self.container.moleculesUnbound()),
			{unboundMolecule}
			)


	@noseAttrib.attr('smalltest', 'transcripts', 'containerObject')
	def test_moleculesBoundOnTranscript(self):
		transcripts = []
		transcriptMolecules = []

		for iTranscript in xrange(2):
			transcript = self.container.transcriptNew(0, 100)
			self.container.transcriptExtend(transcript, 100)

			transcripts.append(transcript)

			positions = [10, 30, 50]
			molecules = []

			for position in positions:
				molecule = self.container.moleculeNew('Ribosome')

				self.container.moleculeLocationIs(molecule, transcript, position, '+',
					10, 5)

				molecules.append(molecule)

			transcriptMolecules.append(molecules)

		for transcript, molecules in zip(transcripts, transcriptMolecules):
			self.assertEqual(
				set(self.container.moleculesBoundOnTranscript(transcript)),
				set(molecules)
				)


	@noseAttrib.attr('smalltest', 'transcripts', 'containerObject')
	def test_moleculesBoundWithName(self):
		transcript = self.container.transcriptNew(0, 100)
		self.container.transcriptExtend(transcript, 100)

		ribosome = self.container.moleculeNew('Ribosome')
		rnase = self.container.moleculeNew('RNAse')

		self.container.moleculeLocationIs(ribosome, transcript, 10, '+',
			10, 5)

		self.container.moleculeLocationIs(rnase, transcript, 30, '+',
			10, 5)

		self.assertEqual(
			set(self.container.moleculesBoundWithName('Ribosome')),
			{ribosome}
			)

		self.assertEqual(
			set(self.container.moleculesBoundWithName('RNAse')),
			{rnase}
			)


def createContainer():
	container = TranscriptsContainer(ARRAY_LENGTH, MOLECULE_ATTRIBUTES)

	return container

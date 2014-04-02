'''

transcripts_container.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 3/19/14

'''


from __future__ import division

import numpy as np

from wholecell.containers.unique_objects_container import UniqueObjectsContainer

MAX_SEARCH_ITERATIONS = 1000 # number of times to search for a place to put a transcript before raising an exception


class TranscriptsContainerException(Exception):
	'''
	TranscriptsContainerException

	An exception subclass raised by the TranscriptsContainer.
	'''
	pass


class TranscriptsContainer(object):
	'''
	TranscriptsContainer

	A container object for transcripts.  Transcripts are contiguous regions of 
	nucleotides, allocated somewhere in a shared array.  Transcripts can be 
	extended, and molecules can be attached to a transcripts.  This 
	implementation is similar to the ChromosomeContainer object, albeit with 
	enough key differences that I've chosen not to inherit from some base 
	implementation.
	'''

	# Special values
	_unused = 0 # a location that is inactive (no binding permitted)
	_reserved = 1 # an inactive location that has been reserved for eventual occupation
	_empty = 2 # an active but empty location

	_inactiveValues = np.array([_unused, _reserved])
	_specialValues = np.array([_unused, _reserved, _empty])
	_offset = _specialValues.size

	# Molecule container properties
	_defaultObjectContainerObjects = {
		'_transcript':{
			'_transPosition':'int64', # location in array
			'_transOrigin':'int64', # origin on chromosome, for determining features/sequence
			'_transExtent':'int64', # length of transcript
			'_transExtentReserved':'int64', # length of region in array reserved
			}
		}

	_defaultObjectContainerAttributes = {
		'_transBound':'bool', # whether or not the molecule is bound to the sequence - after evolveState, all should be bound!
		'_transTranscript':'int64', # ID pointing to the bound transcript
		'_transPosition':'int64', # location in array (not on the transcript)
		'_transDirection':'bool', # False = (+), True = (-)
		'_transExtentForward':'int64', # number of nts
		'_transExtentReverse':'int64', # number of nts
		}

	# Conversion key for directions
	_positiveChar = '+'
	_negativeChar = '-'
	_directionCharToBool = {_positiveChar:False, _negativeChar:True}
	_directionBoolToChar = [_positiveChar, _negativeChar]
	

	def __init__(self, arrayLength, moleculeAttributes, randStream = None):
		self._length = arrayLength

		self._array = np.zeros(self._length, dtype = np.int64)

		molAttrs = self._defaultObjectContainerObjects
		for moleculeName, attributes in moleculeAttributes.viewitems():
			molAttrs[moleculeName] = attributes.copy()
			molAttrs[moleculeName].update(self._defaultObjectContainerAttributes)

		self._objectsContainer =  UniqueObjectsContainer(molAttrs)

		if randStream is None:
			import wholecell.utils.rand_stream

			randStream = wholecell.utils.rand_stream.RandStream()

		self._randStream = randStream


	def transcriptNew(self, chromosomeOrigin, expectedLength = 0):
		position = self._findFreePosition(expectedLength+1)

		transcript = self._objectsContainer.objectNew(
			'_transcript',
			_transPosition = position,
			_transOrigin = chromosomeOrigin,
			_transExtentReserved = expectedLength
			)

		self._array[position] = transcript.attr('_globalIndex') + self._offset

		self._array[position+1:position+1 + expectedLength] = self._reserved

		return transcript


	def transcriptExtend(self, transcript, extent):
		transcriptIndex = transcript.attr('_globalIndex') + self._offset

		position, currentExtent, reserved = transcript.attrs(
			'_transPosition', '_transExtent',
			'_transExtentReserved')

		newExtent = currentExtent + extent

		endPosition = position+1 + newExtent

		region = np.arange(position+1 + currentExtent, endPosition)

		if (newExtent <= reserved) or (
				endPosition < self._length and
				np.setdiff1d(self._array[region], self._inactiveValues).size == 0
				):
			# No risk of collision
			self._array[region] = self._empty

			transcript.attrIs(_transExtent = newExtent)

		else:
			# Find and move to a new location
			oldExtent = np.max([currentExtent, reserved])

			oldRegion = np.arange(position, position+1 + oldExtent)

			oldValues = self._array[oldRegion]
			self._array[oldRegion] = self._unused

			newPosition = self._findFreePosition(np.max([newExtent, oldExtent])+1)

			newRegion = np.arange(newPosition, newPosition+1 + oldExtent)

			self._array[newRegion] = oldValues

			self._array[newPosition+1 + currentExtent:newPosition+1 + newExtent] = self._empty

			transcript.attrIs(
				_transPosition = newPosition,
				_transExtent = newExtent
				)

			positionDelta = newPosition - position

			for molecule in self.moleculesBoundOnTranscript(transcript):
				molecule.attrIs(
					_transPosition = molecule.attr('_transPosition') + positionDelta
					)


	def _findFreePosition(self, extent):
		# TODO: explore better heuristics

		for iteration in xrange(MAX_SEARCH_ITERATIONS):
			position = self._randStream.randi(self._length-extent)

			if (self._array[position:position+extent] == self._unused).all():
				return position

		else:
			raise TranscriptsContainerException(
				'Could not find a free region of size {} in {} iterations'.format(
					extent,
					iteration+1
					)
				)


	def transcriptDel(self, transcript):
		molecules = self.moleculesBoundOnTranscript(transcript)

		for molecule in molecules:
			self.moleculeLocationIsUnbound(molecule)

		position = transcript.attr('_transPosition')
		extent = np.max(transcript.attrs('_transExtent',
				'_transExtentReserved'))

		region = np.arange(position, position+1 + extent)

		self._array[region] = self._unused

		self._objectsContainer.objectDel(transcript)

		return molecules


	def moleculeNew(self, moleculeName, **attributes):
		# Create a new, currently unbound molecule
		return self._objectsContainer.objectNew(moleculeName, **attributes)


	def moleculeDel(self, molecule):
		# Unbind and delete a molecule
		self.moleculeLocationIsUnbound(molecule)
		self._objectsContainer.objectDel(molecule)


	def moleculeLocationIs(self, molecule, transcript, position, direction, extentForward, extentReverse):
		# Set a molecule's location
		transcriptPosition = transcript.attr('_transPosition') + 1
		directionBool = self._directionCharToBool[direction]

		absPosition = position + transcriptPosition

		region = self._region(position + transcriptPosition, directionBool, extentForward, extentReverse)

		moleculeIndex = molecule.attr('_globalIndex') + self._offset

		if np.setdiff1d(self._array[region], [self._empty, moleculeIndex]).size > 0:
			raise TranscriptsContainerException('Attempted to place a molecule in a non-empty region')

		self.moleculeLocationIsUnbound(molecule)

		molecule.attrIs(
			_transBound = True,
			_transTranscript = transcript.attr('_globalIndex') + self._offset,
			_transPosition = absPosition,
			_transDirection = directionBool,
			_transExtentForward = extentForward,
			_transExtentReverse = extentReverse
			)

		self._array[region] = moleculeIndex


	def _region(self, position, directionBool, extentForward, extentReverse):
		# Return a region of a transcript, accounting for directionality
		if directionBool: # == (-)
			return np.arange(position-extentForward+1, position+extentReverse+1)

		else: # == (+)
			return np.arange(position-extentReverse, position+extentForward)


	def moleculeLocation(self, molecule):
		# Return the location of a molecule, if it is bound to a transcript
		if not molecule.attr('_transBound'):
			return None

		else:
			# TODO: implement/decide if this is needed

			raise NotImplementedError()

	# TODO: moleculeTranscript, moleculePosition, moleculeDirection, moleculeFootprint

	def moleculeLocationIsUnbound(self, molecule):
		# Unbind a molecule if it is bound, or do nothing
		if not molecule.attr('_transBound'):
			# Molecule is already unbound
			return

		else:
			# Molecule is bound
			position, directionBool, extentForward, extentReverse = molecule.attrs(
				'_transPosition', '_transDirection', '_transExtentForward', 
				'_transExtentReverse')

			region = self._region(position, directionBool, extentForward, extentReverse)

			self._array[region] = self._empty

			molecule.attrIs(_transBound = False)

	def moleculesBound(self):
		return self._objectsContainer.objects(_transBound = ('==', True))


	def moleculesUnbound(self):
		return self._objectsContainer.objects(_transBound = ('==', False))


	def moleculesBoundOnTranscript(self, transcript):
		transcriptIndex = transcript.attr('_globalIndex') + self._offset
		return self._objectsContainer.objects(_transTranscript = ('==', transcriptIndex))


	def moleculesBoundWithName(self, moleculeName):
		return self._objectsContainer.objectsInCollection(moleculeName, 
			_transBound = ('==', True))


	def moleculeBoundAtPosition(self, transcript, position):
		transcriptPosition = transcript.attr('_transPosition') + 1

		absPosition = position + transcriptPosition

		index = self._array[absPosition] - self._offset

		if index < 0:
			return None

		else:
			return self._objectsContainer._objectByGlobalIndex(index)


	def moleculesBoundOverExtent(self, transcript, position, direction, extentForward, extentReverse):
		transcriptPosition = transcript.attr('_transPosition') + 1

		absPosition = position + transcriptPosition

		directionBool = self._directionCharToBool[direction]

		region = self._region(absPosition, directionBool, extentForward, extentReverse)

		indexes = np.setdiff1d(self._array[region], self._specialValues) - self._offset

		return self._objectsContainer._objectsByGlobalIndex(indexes)


	def moleculesBoundNearMolecule(self, molecule, extentForward, extentReverse):
		raise NotImplementedError()


	def __eq__(self, other):
		return (self._array == other._array).all()

	# TODO: save/load

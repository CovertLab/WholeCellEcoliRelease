'''

transcripts_container.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 3/19/14

'''


from __future__ import division

import numpy as np

from wholecell.containers.unique_objects_container import UniqueObjectsContainer


class TranscriptsContainerException(Exception):
	'''
	TranscriptsContainerException

	An exception subclass raised by the TranscriptsContainer.
	'''
	pass


class TranscriptsContainer(object):
	'''
	TranscriptsContainer

	TODO: doc
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
			'_transcriptPosition':'int64', # location in array
			'_transcriptOrigin':'int64', # origin on chromosome, for determining features/sequence
			'_transcriptExtent':'int64', # length of transcript
			'_transcriptExtentReserved':'int64', # length of region in array reserved
			}
		}

	_defaultObjectContainerAttributes = {
		'_transcriptBound':'bool', # whether or not the molecule is bound to the sequence - after evolveState, all should be bound!
		'_transcript':'int64', # ID pointing to the bound transcript
		'_transcriptPosition':'int64', # location in array (not on the transcript)
		'_transcriptDirection':'bool', # False = (+), True = (-)
		'_transcriptExtentForward':'int64', # number of nts
		'_transcriptExtentReverse':'int64', # number of nts
		}

	# Conversion key for directions
	_positiveChar = '+'
	_negativeChar = '-'
	_directionCharToBool = {_positiveChar:False, _negativeChar:True}
	_directionBoolToChar = [_positiveChar, _negativeChar]
	

	def __init__(self, arrayLength):
		self._length = arrayLength

		self._array = np.zeros(self._length, dtype = np.int64)

		molAttrs = self._defaultObjectContainerObjects
		for moleculeName, attributes in moleculeAttributes.viewitems():
			molAttrs[moleculeName] = attributes.copy()
			molAttrs[moleculeName].update(self._defaultObjectContainerAttributes)

		self._objectsContainer =  UniqueObjectsContainer(molAttrs)


	def transcriptNew(self, chromosomeOrigin, expectedLength = 0):
		position = self._findFreePosition(expectedLength+1)

		transcript = self._objectsContainer.objectNew(
			'_transcript',
			_transcriptPosition = position,
			_transcriptOrigin = chromosomeOrigin,
			_transcriptExtentReserved = expectedLength
			)

		self._array[position] = transcript.attr('_globalIndex') + self._offset

		self._array[position+1:position+1+expectedLength] = self._reserved


	def transcriptExtend(self, transcript, extent):
		position = transcript.attr('_transcriptPosition')
		transcriptIndex = transcript.attr('_globalIndex') + self._offset
		currentLength = transcript.attr('_transcriptExtent')
		reserved = transcript.attr('_transcriptExtentReserved')

		newLength = currentLength + extent

		endPosition = position+1 + newLength

		region = np.arange(position+1 + currentLength, endPosition)

		if (newLength <= reserved) or (
				np.setdiff1d(self._array[region], self._inactiveValues).size == 0
				and endPosition < self._length
				):
			# No risk of collision
			self._array[region] = self._empty

		else:
			# Find and move to a new location
			oldRegion = np.arange(position, position+1+currentLength)

			oldValues = self._array[oldRegion]
			self._array[oldRegion] = self._unused

			newPosition = self._findFreePosition(newLength+1)

			# TODO: set new region to old values
			# TODO: update transcript and molecules with new positions

			raise NotImplementedError()


	def _findFreePosition(self, extent):
		raise NotImplementedError()


	def transcriptDel(self, transcript):
		raise NotImplementedError()


	def moleculeNew(self, moleculeName, **attributes):
		# Create a new, currently unbound molecule
		return self._objectsContainer.objectNew(moleculeName, **attributes)


	def moleculeDel(self, molecule):
		# Unbind and delete a molecule
		self.moleculeLocationIsUnbound(molecule)
		self._objectsContainer.objectDel(molecule)


	def moleculeLocationIs(self, molecule, transcript, position, direction, extentForward, extentReverse):
		# Set a molecule's location
		strandIndex = self._strandNameToIndex[strand]
		directionBool = self._directionCharToBool[direction]

		region = self._region(position, directionBool, extentForward, extentReverse)

		if not (self._array[strandIndex, region] == self._empty).all():
			raise TranscriptsContainerException('Attempted to place a molecule in a non-empty region')

		self.moleculeLocationIsUnbound(molecule)

		molecule.attrIs('_transcriptBound', True) # TODO: attrsAre method
		molecule.attrIs('_transcriptStrand', strandIndex)
		molecule.attrIs('_transcriptPosition', position)
		molecule.attrIs('_transcriptDirection', directionBool)
		molecule.attrIs('_transcriptExtentForward', extentForward)
		molecule.attrIs('_transcriptExtentReverse', extentReverse)

		self._array[strandIndex, region] = molecule.attr('_globalIndex') + self._offset


	def _region(self, position, directionBool, extentForward, extentReverse):
		# Return a region of a strand, accounting for circularity and directionality
		if directionBool: # == (-)
			return np.arange(position-extentForward+1, position+extentReverse+1) % self._length

		else: # == (+)
			return np.arange(position-extentReverse, position+extentForward) % self._length


	def moleculeLocation(self, molecule):
		# Return the location a molecule, if it is bound to the chromosome
		if not molecule.attr('_transcriptBound'):
			return None

		else:
			return (
				self._strandNames[molecule.attr('_transcriptStrand')], # TODO: attrs method
				molecule.attr('_transcriptPosition'),
				self._directionBoolToChar[molecule.attr('_transcriptDirection')],
				molecule.attr('_transcriptExtentForward'),
				molecule.attr('_transcriptExtentReverse'),
				)

	# TODO: moleculeTranscript, moleculePosition, moleculeDirection, moleculeFootprint

	def moleculeLocationIsUnbound(self, molecule):
		# Unbind a molecule (from a normal location or a fork) if it is bound, or do nothing
		if not molecule.attr('_transcriptBound'):
			# Molecule is already unbound
			return

		else:
			# Molecule is bound to a normal location (not a fork)
			strandIndex = molecule.attr('_transcriptStrand')
			position = molecule.attr('_transcriptPosition')
			directionBool = molecule.attr('_transcriptDirection')
			extentForward = molecule.attr('_transcriptExtentForward')
			extentReverse = molecule.attr('_transcriptExtentReverse')

			region = self._region(position, directionBool, extentForward, extentReverse)

			self._array[strandIndex, region] = self._empty

			molecule.attrIs('_transcriptBound', False)


	def moleculesBound(self, moleculeName = None, transcript = None, 
			position = None, direction = None, extentForward = None,
			extentReverse = None):
		# Returns bound molecules, with sets of optional arguments:
		# moleculeName: only molecules with this name
		# transcript: molecules on this transcript
		# +position, direction: molecule at a specific position
		# +extentForward, extentReverse: molecules over a region

		# TODO: check for inconsistent sets of arguments
		# TODO: make the queries more efficient

		specifiesName = moleculeName is not None

		specifiesPosition = (strand is not None and position is not None)

		specifiesExtent = (direction is not None and extentForward is not None
			and extentReverse is not None)

		specificRequest = specifiesName or specifiesPosition

		if specifiesPosition:
			strandIndex = self._strandNameToIndex[strand]

			if specifiesExtent:
				directionBool = self._directionCharToBool[direction]

				region = self._region(position, directionBool, extentForward, extentReverse)

				indexes = np.setdiff1d(
					self._array[strandIndex, region],
					self._specialValues
					) - self._offset

			else:
				indexes = np.setdiff1d(
					self._array[strandIndex, position],
					self._specialValues
					) - self._offset

		else:
			indexes = np.setdiff1d(self._array, self._specialValues) - self._offset

		if specifiesName:
			molecules = self._objectsContainer._objectsByGlobalIndex(indexes)

			return {molecule for molecule in molecules
				if molecule.name() == moleculeName}

		else:
			return self._objectsContainer._objectsByGlobalIndex(indexes)


	def moleculesUnbound(self):
		# Returns a set of all unbound molecules
		raise NotImplementedError()


	def __eq__(self, other):
		return (self._array == other._array).all()

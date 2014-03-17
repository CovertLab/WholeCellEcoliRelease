
from __future__ import division

import numpy as np

import wholecell.utils.unique_objects_container

class ChromosomeBoundMoleculeContainer(object):
	# Special values for the state of the chromosome
	_inactive = 0 # a location that is inactive (no binding permitted)
	_empty = 1 # an active but empty location

	_specialValues = np.array([_inactive, _empty])
	_offset = _specialValues.size

	# Molecule container properties
	_defaultObjectContainerObjects = {
		'_fork':{
			'_sequenceStrand':'uint32', # index of strand
			'_sequencePosition':'uint32', # bound nt
			'_sequenceDirection':'bool', # False = (+), True = (-)
			}
		}

	_defaultObjectContainerAttributes = {
		'_sequenceBound':'bool', # whether or not the molecule is bound to the sequence - after evolveState, all should be bound!
		'_sequenceStrand':'uint32', # index of strand
		'_sequencePosition':'uint32', # bound nt
		'_sequenceDirection':'bool', # False = (+), True = (-)
		'_sequenceExtentForward':'uint32', # number of nts
		'_sequenceExtentReverse':'uint32', # number of nts
		'_sequenceBoundToFork':'bool', # special property for molecules bound on a forked region
		}

	# Single-character names for denoting strand identity
	_rootChar = 'R' 
	_childAChar = 'A'
	_childBChar = 'B'

	# Conversion key for directions
	_positiveChar = '+'
	_negativeChar = '-'
	_directionCharToBool = {_positiveChar:False, _negativeChar:True}
	_directionBoolToChar = [_positiveChar, _negativeChar]
	
	def __init__(self, nBases, strandMultiplicity, moleculeAttributes):
		self._length = nBases

		self._strandMultiplicity = strandMultiplicity
		self._buildStrandConnectivity()
		
		self._array = np.zeros((self._nStrands, self._length), dtype = np.int32) # TODO: choose best dtype based on array size

		self._array[0, :] = self._empty # Root strand is always active

		molAttrs = self._defaultObjectContainerObjects
		for moleculeName, attributes in moleculeAttributes.viewitems():
			molAttrs[moleculeName] = attributes.copy()
			molAttrs[moleculeName].update(self._defaultObjectContainerAttributes)

		self._objectsContainer =  wholecell.utils.unique_objects_container.UniqueObjectsContainer(
			molAttrs)


	def _buildStrandConnectivity(self):
		self._strandNames = []
		self._strandNames.append(self._rootChar)

		parentNames = self._strandNames
		for i in xrange(self._strandMultiplicity - 1):
			childNames = []

			for parentName in parentNames:
				childNames.append(parentName + self._childAChar)
				childNames.append(parentName + self._childBChar)

			self._strandNames.extend(childNames)
			parentNames = childNames
		
		self._nStrands = len(self._strandNames)

		self._strandNameToIndex = {name:ind for ind, name in enumerate(self._strandNames)}
		self._strandChildrenIndexes = [
			(self._strandNameToIndex[name + self._childAChar], self._strandNameToIndex[name + self._childBChar])
			if self._strandNameToIndex.has_key(name + self._childAChar) else None
			for name in self._strandNames
			]

		self._strandParentIndex = [None]*self._nStrands

		for parentIndex, childrenIndexes in enumerate(self._strandChildrenIndexes):
			if childrenIndexes is not None:
				self._strandParentIndex[childrenIndexes[0]] = parentIndex
				self._strandParentIndex[childrenIndexes[1]] = parentIndex


	def moleculeNew(self, moleculeName, **attributes):
		return self._objectsContainer.objectNew(moleculeName, **attributes)


	def moleculeDel(self, molecule):
		self.moleculeLocationIsUnbound(molecule)
		self._objectsContainer.objectDel(molecule)


	def moleculeLocationIs(self, molecule, strand, position, direction, extentForward, extentReverse):
		strandIndex = self._strandNameToIndex[strand]
		directionBool = self._directionCharToBool[direction]

		extentPositive, extentNegative = self._extentRelativeToAbsolute(
			extentForward, extentReverse, directionBool)

		region = self._region(strandIndex, position, directionBool, extentPositive, extentNegative)

		assert (region == self._empty).all(), 'Attempted to place a molecule in a non-empty region'

		self.moleculeLocationIsUnbound(molecule)

		molecule.attrIs('_sequenceBound', True) # TODO: attrsAre method
		molecule.attrIs('_sequenceStrand', strandIndex)
		molecule.attrIs('_sequencePosition', position)
		molecule.attrIs('_sequenceDirection', directionBool)
		molecule.attrIs('_sequenceExtentForward', extentForward)
		molecule.attrIs('_sequenceExtentReverse', extentReverse)

		region[:] = molecule.attr('_globalIndex') + self._offset


	def moleculeLocationIsFork(self, molecule, fork, extentForward, extentReverse): # TODO: different child extents
		# NOTE: molecule orientation assumed to be in the direction of the fork

		assert (extentForward > 0 or extentReverse > 0), 'The footprint of a molecule placed on a fork must be > 1'

		forkStrand = fork.attr('_sequenceStrand')
		forkPosition = fork.attr('_sequencePosition')
		forkDirection = fork.attr('_sequenceDirection')

		(regionParent, regionChildA, regionChildB) = self._forkedRegions(
			forkStrand, forkPosition, forkDirection, extentForward,
			extentReverse) 

		assert (regionParent == self._empty).all(), 'Attempted to place a molecule in a non-empty region'
		assert (regionChildA == self._empty).all(), 'Attempted to place a molecule in a non-empty region'
		assert (regionChildB == self._empty).all(), 'Attempted to place a molecule in a non-empty region'

		self.moleculeLocationIsUnbound(molecule)

		molecule.attrIs('_sequenceBound', True)
		molecule.attrIs('_sequenceStrand', forkStrand)
		molecule.attrIs('_sequencePosition', forkPosition)
		molecule.attrIs('_sequenceDirection', forkDirection)
		molecule.attrIs('_sequenceExtentForward', extentForward)
		molecule.attrIs('_sequenceExtentReverse', extentReverse)
		molecule.attrIs('_sequenceBoundToFork', True)

		index = molecule.attr('_globalIndex') + self._offset

		regionParent[:] = index
		regionChildA[:] = index
		regionChildB[:] = index


	def _region(self, strandIndex, position, directionBool, extentForward, extentReverse):
		extentPositive, extentNegative = self._extentRelativeToAbsolute(
			extentForward, extentReverse, directionBool)

		return self._array[
			strandIndex,
			np.arange(position-extentNegative, position+extentPositive) % self._length
			]		


	def _forkedRegions(self, forkStrand, forkPosition, forkDirection, extentForward, extentReverse):
		(childStrandA, childStrandB) = self._strandChildrenIndexes[forkStrand]

		extentPositive, extentNegative = self._extentRelativeToAbsolute(
			extentForward, extentReverse, forkDirection)

		if forkDirection: # True = (-)
			regionParent = self._array[forkStrand, forkPosition-extentNegative+1:forkPosition]
			regionChildA = self._array[childStrandA, forkPosition:forkPosition+extentPositive]
			regionChildB = self._array[childStrandB, forkPosition:forkPosition+extentPositive]

		else:
			regionParent = self._array[forkStrand, forkPosition+1:forkPosition+extentPositive]
			regionChildA = self._array[childStrandA, forkPosition-extentNegative+1:forkPosition+1]
			regionChildB = self._array[childStrandB, forkPosition-extentNegative+1:forkPosition+1]

		return (regionParent, regionChildA, regionChildB)


	# TODO: assertions about correct footprints after operations

	def moleculeLocation(self, molecule):
		if not molecule.attr('_sequenceBound'):
			return None

		else:
			return (
				self._strandNames[molecule.attr('_sequenceStrand')], # TODO: attrs method
				molecule.attr('_sequencePosition'),
				self._directionBoolToChar[molecule.attr('_sequenceDirection')],
				molecule.attr('_sequenceExtentForward'),
				molecule.attr('_sequenceExtentReverse'),
				)

	# TODO: moleculeStrand, moleculePosition, moleculeDirection, moleculeFootprint

	def moleculeLocationIsUnbound(self, molecule):
		if not molecule.attr('_sequenceBound'):
			return


		elif molecule.attr('_sequenceBoundToFork'):
			strandIndex = molecule.attr('_sequenceStrand')
			position = molecule.attr('_sequencePosition')
			directionBool = molecule.attr('_sequenceDirection')
			extentForward = molecule.attr('_sequenceExtentForward')
			extentReverse = molecule.attr('_sequenceExtentReverse')

			(regionParent, regionChildA, regionChildB) = self._forkedRegions(
				strandIndex, position, directionBool, extentForward, extentReverse)

			regionParent[:] = self._empty
			regionChildA[:] = self._empty
			regionChildB[:] = self._empty

			molecule.attrIs('_sequenceBound', False)
			molecule.attrIs('_sequenceBoundToFork', False)


		else:
			strandIndex = molecule.attr('_sequenceStrand')
			position = molecule.attr('_sequencePosition')
			directionBool = molecule.attr('_sequenceDirection')
			extentForward = molecule.attr('_sequenceExtentForward')
			extentReverse = molecule.attr('_sequenceExtentReverse')

			extentPositive, extentNegative = self._extentRelativeToAbsolute(
				extentForward, extentReverse, directionBool)

			region = self._region(strandIndex, position, directionBool, extentPositive, extentNegative)

			region[:] = self._empty

			molecule.attrIs('_sequenceBound', False)


	def _extentRelativeToAbsolute(self, extentForward, extentReverse, directionBool):
		if directionBool: # flip if negative
			return extentReverse, extentForward

		else:
			return extentForward, extentReverse


	def moleculesBound(self, moleculeName = None, strand = None, 
			position = None, direction = None, extentForward = None,
			extentReverse = None):

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

				extentPositive, extentNegative = self._extentRelativeToAbsolute(
					extentForward, extentReverse, directionBool)

				region = self._region(strandIndex, position, directionBool, extentPositive, extentNegative)

				indexes = np.setdiff1d(region) - self._offset

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
		raise NotImplementedError()


	def divideRegion(self, strandName, start, stop): # TODO: handle stop < start
		# NOTE: start to stop is inclusive, unlike "range"!
		strandParent = self._strandNameToIndex[strandName]
		try:
			strandChildA, strandChildB = self._strandChildrenIndexes[strandParent]

		except TypeError:
			raise Exception('No space allocated for strand {} to divide into'.format(strandName))


		assert (self._array[strandParent, start:stop+1] == self._empty).all(), 'Attempted to divide a non-empty or non-existent region'

		self._array[strandParent, start:stop+1] = self._inactive

		self._array[strandChildA, start:stop+1] = self._empty
		self._array[strandChildB, start:stop+1] = self._empty

		forkStart = self._objectsContainer.objectNew(
			'_fork',
			_sequenceStrand = strandParent,
			_sequencePosition = start,
			_sequenceDirection = 1, # NOTE: I've chosen the convention that the fork "direction" is towards the parent (i.e. the natural direction of extension)
			)

		forkStop = self._objectsContainer.objectNew(
			'_fork',
			_sequenceStrand = strandParent,
			_sequencePosition = stop,
			_sequenceDirection = 0,
			)

		self._array[strandParent, start] = forkStart.attr('_globalIndex') + self._offset
		self._array[strandParent, stop] = forkStop.attr('_globalIndex') + self._offset

		return forkStart, forkStop


	def forks(self, strand = None, position = None, direction = None,
			extentForward = None, extentReverse = None):

		return self.moleculesBound('_fork', strand, position, direction,
			extentForward, extentReverse)


	def forkExtend(self, fork, extent):
		forkStrand = fork.attr('_sequenceStrand')
		forkPosition = fork.attr('_sequencePosition')
		forkDirection = fork.attr('_sequenceDirection')

		strandChildA, strandChildB = self._strandChildrenIndexes[forkStrand]

		if forkDirection == 0:
			region = np.arange(forkPosition+1, forkPosition+extent+1) % self._length
			newPosition = (forkPosition+extent) % self._length

		else:
			region = np.arange(forkPosition-extent, forkPosition) % self._length
			newPosition = (forkPosition-extent) % self._length

		assert (self._array[forkStrand, region] == self._empty).all(), 'Attempted to extend a fork into a non-empty region'

		self._array[forkStrand, region] = self._inactive

		self._array[strandChildA, region] = self._empty
		self._array[strandChildB, region] = self._empty

		self._array[forkStrand, forkPosition] = self._inactive
		fork.attrIs('_sequencePosition', newPosition)
		self._array[forkStrand, newPosition] = fork.attr('_globalIndex') + self._offset

		return newPosition


	def forksCombine(self, fork1, fork2):
		# Combine two forks on the same strand, splitting the chromosome
		raise NotImplementedError()


	def moleculeOnFork(self, fork):
		raise NotImplementedError()


	# def findLocationToBind(self, width):
	# 	return self.findLocationsToBind(width, 1)[0]


	# def findLocationsToBind(self, width, nLocations):
	# 	# TODO: make this method faster (cache valid locations by footprint, use sparse representation...)

	# 	# locations = 
	# 	pass


	# TODO: use this method
	# def _range(self, start, stop):
	# 	# Converts a start:stop slice into a range that handles the circularity
	# 	# of the chromosome

	# 	return np.arange(start, stop) % self._length
	

	# TODO: circularly-permuted indexing
	# TODO: saving
	# TODO: update container time, flush deleted molecules, update queries?
	# TODO: handle/pass sequence, multiplicity
	# TODO: write tests
'''

chromosome_container.py

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@data: Created 3/12/14

'''


from __future__ import division

import numpy as np

from wholecell.containers.unique_objects_container import UniqueObjectsContainer


class ChrosomeContainerException(Exception):
	'''
	ChrosomeContainerException

	An exception subclass raised by the ChromosomeContainer.
	'''
	pass


class ChromosomeContainer(object):
	'''
	ChromosomeContainer

	A container class.  Tracks the locations of objects (molecules) on a
	circular, branched, directional structure, with provided footprints.  
	Objects can be attached to forks as well, and forks can be extended along 
	their parent strand.
	'''

	# Special values for the state of the chromosome
	_inactive = 0 # a location that is inactive (no binding permitted)
	_empty = 1 # an active but empty location

	_specialValues = np.array([_inactive, _empty])
	_offset = _specialValues.size

	# Molecule container properties
	_defaultObjectContainerObjects = {
		'_fork':{
			'_sequenceStrand':'int64', # index of strand
			'_sequencePosition':'int64', # bound nt
			'_sequenceDirection':'bool', # False = (+), True = (-)
			}
		}

	_defaultObjectContainerAttributes = {
		'_sequenceBound':'bool', # whether or not the molecule is bound to the sequence - after evolveState, all should be bound!
		'_sequenceStrand':'int64', # index of strand
		'_sequencePosition':'int64', # bound nt
		'_sequenceDirection':'bool', # False = (+), True = (-)
		'_sequenceExtentForward':'int64', # number of nts
		'_sequenceExtentReverse':'int64', # number of nts
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
		
		self._array = np.zeros((self._nStrands, self._length), dtype = np.int64)

		self._array[0, :] = self._empty # Root strand is always active

		molAttrs = self._defaultObjectContainerObjects
		for moleculeName, attributes in moleculeAttributes.viewitems():
			molAttrs[moleculeName] = attributes.copy()
			molAttrs[moleculeName].update(self._defaultObjectContainerAttributes)

		self._objectsContainer =  UniqueObjectsContainer(molAttrs)


	def _buildStrandConnectivity(self):
		# Build the lists and tables needed for interactions between strands

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
		# Create a new, currently unbound molecule
		return self._objectsContainer.objectNew(moleculeName, **attributes)


	def moleculeDel(self, molecule):
		# Unbind and delete a molecule
		self.moleculeLocationIsUnbound(molecule)
		self._objectsContainer.objectDel(molecule)


	def moleculeLocationIs(self, molecule, strand, position, direction, extentForward, extentReverse):
		# Set a molecule's location
		strandIndex = self._strandNameToIndex[strand]
		directionBool = self._directionCharToBool[direction]

		region = self._region(position, directionBool, extentForward, extentReverse)

		if not (self._array[strandIndex, region] == self._empty).all():
			raise ChrosomeContainerException('Attempted to place a molecule in a non-empty region')

		self.moleculeLocationIsUnbound(molecule)

		molecule.attrIs('_sequenceBound', True) # TODO: attrsAre method
		molecule.attrIs('_sequenceStrand', strandIndex)
		molecule.attrIs('_sequencePosition', position)
		molecule.attrIs('_sequenceDirection', directionBool)
		molecule.attrIs('_sequenceExtentForward', extentForward)
		molecule.attrIs('_sequenceExtentReverse', extentReverse)

		self._array[strandIndex, region] = molecule.attr('_globalIndex') + self._offset


	def moleculeLocationIsFork(self, molecule, fork, extentForward, extentReverse): # TODO: different child extents
		# Assign a molecule's location to that of a fork
		# NOTE: molecule orientation assumed to be in the direction of the fork

		if not (extentForward > 0 or extentReverse > 0):
			raise ChrosomeContainerException('The footprint of a molecule placed on a fork must be > 1')

		forkStrand = fork.attr('_sequenceStrand')
		forkPosition = fork.attr('_sequencePosition')
		forkDirection = fork.attr('_sequenceDirection')

		(regionParent, regionChildA, regionChildB) = self._forkedRegions(
			forkPosition, forkDirection, extentForward,	extentReverse)

		(childStrandA, childStrandB) = self._strandChildrenIndexes[forkStrand]

		if not (self._array[forkStrand, regionParent] == self._empty).all():
			raise ChrosomeContainerException('Attempted to place a molecule in a non-empty region')

		if not (self._array[childStrandA, regionChildA] == self._empty).all():
			raise ChrosomeContainerException('Attempted to place a molecule in a non-empty region')

		if not (self._array[childStrandB, regionChildB] == self._empty).all():
			raise ChrosomeContainerException('Attempted to place a molecule in a non-empty region')

		self.moleculeLocationIsUnbound(molecule)

		molecule.attrIs('_sequenceBound', True)
		molecule.attrIs('_sequenceStrand', forkStrand)
		molecule.attrIs('_sequencePosition', forkPosition)
		molecule.attrIs('_sequenceDirection', forkDirection)
		molecule.attrIs('_sequenceExtentForward', extentForward)
		molecule.attrIs('_sequenceExtentReverse', extentReverse)
		molecule.attrIs('_sequenceBoundToFork', True)

		index = molecule.attr('_globalIndex') + self._offset

		self._array[forkStrand, regionParent] = index
		self._array[childStrandA, regionChildA] = index
		self._array[childStrandB, regionChildB] = index


	def _region(self, position, directionBool, extentForward, extentReverse):
		# Return a region of a strand, accounting for circularity and directionality
		if directionBool: # == (-)
			return np.arange(position-extentForward+1, position+extentReverse+1) % self._length

		else: # == (+)
			return np.arange(position-extentReverse, position+extentForward) % self._length


	def _forkedRegions(self, forkPosition, forkDirection, extentForward, extentReverse):
		# Return regions on strands surrounding a fork, accounting for circularity and directionality
		if forkDirection: # == (-)
			regionParent = np.arange(forkPosition-extentForward+1, forkPosition) % self._length
			regionChildA = np.arange(forkPosition, forkPosition+extentReverse+1) % self._length
			regionChildB = np.arange(forkPosition, forkPosition+extentReverse+1) % self._length

		else: # == (+)
			regionParent = np.arange(forkPosition+1, forkPosition+extentForward) % self._length
			regionChildA = np.arange(forkPosition-extentReverse, forkPosition+1) % self._length
			regionChildB = np.arange(forkPosition-extentReverse, forkPosition+1) % self._length

		return (regionParent, regionChildA, regionChildB)


	def moleculeLocation(self, molecule):
		# Return the location a molecule, if it is bound to the chromosome
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
		# Unbind a molecule (from a normal location or a fork) if it is bound, or do nothing
		if not molecule.attr('_sequenceBound'):
			# Molecule is already unbound
			return


		elif molecule.attr('_sequenceBoundToFork'):
			# Molecule is bound to a fork
			strandIndex = molecule.attr('_sequenceStrand')
			position = molecule.attr('_sequencePosition')
			directionBool = molecule.attr('_sequenceDirection')
			extentForward = molecule.attr('_sequenceExtentForward')
			extentReverse = molecule.attr('_sequenceExtentReverse')

			(childStrandA, childStrandB) = self._strandChildrenIndexes[strandIndex]

			(regionParent, regionChildA, regionChildB) = self._forkedRegions(
				position, directionBool, extentForward, extentReverse)

			self._array[strandIndex, regionParent] = self._empty
			self._array[childStrandA, regionChildA] = self._empty
			self._array[childStrandB, regionChildB] = self._empty

			molecule.attrIs('_sequenceBound', False)
			molecule.attrIs('_sequenceBoundToFork', False)


		else:
			# Molecule is bound to a normal location (not a fork)
			strandIndex = molecule.attr('_sequenceStrand')
			position = molecule.attr('_sequencePosition')
			directionBool = molecule.attr('_sequenceDirection')
			extentForward = molecule.attr('_sequenceExtentForward')
			extentReverse = molecule.attr('_sequenceExtentReverse')

			region = self._region(position, directionBool, extentForward, extentReverse)

			self._array[strandIndex, region] = self._empty

			molecule.attrIs('_sequenceBound', False)


	def moleculesBound(self, moleculeName = None, strand = None, 
			position = None, direction = None, extentForward = None,
			extentReverse = None):
		# Returns bound molecules, with sets of optional arguments:
		# moleculeName: only molecules with this name
		# strand, position, direction: molecule at a specific position
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


	def divideRegion(self, strandName, start, stop):
		# Break an empty, continuous region of a strand into two child strands,
		# and return the two "fork" objects

		# NOTE: start to stop is inclusive, unlike "range"!
		strandParent = self._strandNameToIndex[strandName]
		try:
			strandChildA, strandChildB = self._strandChildrenIndexes[strandParent]

		except TypeError:
			raise ChrosomeContainerException('No space allocated for strand {} to divide into'.format(strandName))

		# raise exception if start/stop outside length
		# raise exception if start == stop

		if stop < start:
			region = np.r_[np.arange(start, self._length), np.arange(stop+1)]

		else:
			region = np.arange(start, stop+1)

		if not (self._array[strandParent, region] == self._empty).all():
			raise ChrosomeContainerException('Attempted to divide a non-empty or non-existent region')

		self._array[strandParent, region] = self._inactive

		self._array[strandChildA, region] = self._empty
		self._array[strandChildB, region] = self._empty

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
		# Return a set of forks

		return self.moleculesBound('_fork', strand, position, direction,
			extentForward, extentReverse)


	def forkExtend(self, fork, extent):
		# Move a fork along its parent strand (shortening the parent and 
		# elongating the children), and return the new position

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

		if not (self._array[forkStrand, region] == self._empty).all():
			raise ChrosomeContainerException('Attempted to extend a fork into a non-empty region')

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
		# Return the molecule on a fork, if any
		raise NotImplementedError()


	def rootStrand(self):
		# Return the character used to identify the root strand
		return self._rootChar


	#def childStrands
	#def parentStrand


	def __eq__(self, other):
		return (self._array == other._array).all()


	# def findLocationToBind(self, width):
	# 	return self.findLocationsToBind(width, 1)[0]


	# def findLocationsToBind(self, width, nLocations):
	# 	# TODO: make this method faster (cache valid locations by footprint, use sparse representation...)

	# 	# locations = 
	# 	pass
	

	# TODO: saving
	# TODO: update container time, flush deleted molecules, update queries?
	# TODO: handle/pass sequence

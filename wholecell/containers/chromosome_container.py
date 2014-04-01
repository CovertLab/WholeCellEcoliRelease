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
			'_chromStrand':'int64', # index of strand
			'_chromPosition':'int64', # bound nt
			'_chromDirection':'bool', # False = (+), True = (-)
			}
		}

	_defaultObjectContainerAttributes = {
		'_chromBound':'bool', # whether or not the molecule is bound to the sequence - after evolveState, all should be bound!
		'_chromStrand':'int64', # index of strand
		'_chromPosition':'int64', # bound nt
		'_chromDirection':'bool', # False = (+), True = (-)
		'_chromExtentForward':'int64', # number of nts
		'_chromExtentReverse':'int64', # number of nts
		'_chromBoundToFork':'bool', # special property for molecules bound on a forked region
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

		moleculeIndex = molecule.attr('_globalIndex') + self._offset

		if np.setdiff1d(self._array[strandIndex, region], [self._empty, moleculeIndex]).size > 0:
			raise ChrosomeContainerException('Attempted to place a molecule in a non-empty region')

		self.moleculeLocationIsUnbound(molecule)

		molecule.attrIs(
			_chromBound = True,
			_chromStrand = strandIndex,
			_chromPosition = position,
			_chromDirection = directionBool,
			_chromExtentForward = extentForward,
			_chromExtentReverse = extentReverse
			)

		self._array[strandIndex, region] = moleculeIndex


	def moleculeLocationIsFork(self, molecule, fork, extentForward, extentReverse): # TODO: different child extents
		# Assign a molecule's location to that of a fork
		# NOTE: molecule orientation assumed to be in the direction of the fork

		if not (extentForward > 0 or extentReverse > 0):
			raise ChrosomeContainerException('The footprint of a molecule placed on a fork must be > 1')

		forkStrand, forkPosition, forkDirection = fork.attrs('_chromStrand',
			'_chromPosition', '_chromDirection')

		(regionParent, regionChildA, regionChildB) = self._forkedRegions(
			forkPosition, forkDirection, extentForward,	extentReverse)

		(childStrandA, childStrandB) = self._strandChildrenIndexes[forkStrand]

		moleculeIndex = molecule.attr('_globalIndex') + self._offset

		if np.setdiff1d(self._array[forkStrand, regionParent], [self._empty, moleculeIndex]).size > 0:
			raise ChrosomeContainerException('Attempted to place a molecule in a non-empty region')

		if np.setdiff1d(self._array[childStrandA, regionChildA], [self._empty, moleculeIndex]).size > 0:
			raise ChrosomeContainerException('Attempted to place a molecule in a non-empty region')

		if np.setdiff1d(self._array[childStrandB, regionChildB], [self._empty, moleculeIndex]).size > 0:
			raise ChrosomeContainerException('Attempted to place a molecule in a non-empty region')

		self.moleculeLocationIsUnbound(molecule)

		molecule.attrIs(
			_chromBound = True,
			_chromStrand = forkStrand,
			_chromPosition = forkPosition,
			_chromDirection = forkDirection,
			_chromExtentForward = extentForward,
			_chromExtentReverse = extentReverse,
			_chromBoundToFork = True
			)

		index = molecule.attr('_globalIndex') + self._offset

		self._array[forkStrand, regionParent] = index
		self._array[childStrandA, regionChildA] = index
		self._array[childStrandB, regionChildB] = index


	def _region(self, position, directionBool, extentForward, extentReverse):
		# Return a region of a strand, accounting for circularity and directionality
		if directionBool: # == (-)
			return np.arange(position+extentReverse, position-extentForward, -1) % self._length

		else: # == (+)
			return np.arange(position-extentReverse, position+extentForward) % self._length


	def _forkedRegions(self, forkPosition, forkDirection, extentForward, extentReverse):
		# Return regions on strands surrounding a fork, accounting for circularity and directionality
		if forkDirection: # == (-)
			# regionParent = np.arange(forkPosition-extentForward+1, forkPosition) % self._length

			regionParent = np.arange(forkPosition-1, forkPosition-extentForward, -1) % self._length
			regionChildA = np.arange(forkPosition, forkPosition+extentReverse+1) % self._length
			regionChildB = np.arange(forkPosition, forkPosition+extentReverse+1) % self._length

		else: # == (+)
			regionParent = np.arange(forkPosition+1, forkPosition+extentForward) % self._length

			# regionChildA = np.arange(forkPosition-extentReverse, forkPosition+1) % self._length
			# regionChildB = np.arange(forkPosition-extentReverse, forkPosition+1) % self._length

			regionChildA = np.arange(forkPosition, forkPosition-extentReverse-1, -1) % self._length
			regionChildB = np.arange(forkPosition, forkPosition-extentReverse-1, -1) % self._length

		return (regionParent, regionChildA, regionChildB)


	def moleculeLocation(self, molecule):
		# Return the location a molecule, if it is bound to the chromosome
		if not molecule.attr('_chromBound'):
			return None

		else:
			return (
				self._strandNames[molecule.attr('_chromStrand')], # TODO: attrs method
				molecule.attr('_chromPosition'),
				self._directionBoolToChar[molecule.attr('_chromDirection')],
				molecule.attr('_chromExtentForward'),
				molecule.attr('_chromExtentReverse'),
				)

	# TODO: moleculeStrand, moleculePosition, moleculeDirection, moleculeFootprint


	def moleculeFootprint(self, molecule):
		position, direction, extentForward, extentReverse = molecule.attrs(
			'_chromPosition', '_chromDirection', '_chromExtentForward',
			'_chromExtentReverse')

		return self._region(position, direction, extentForward, extentReverse)


	def moleculeLocationIsUnbound(self, molecule):
		# Unbind a molecule (from a normal location or a fork) if it is bound, or do nothing
		if not molecule.attr('_chromBound'):
			# Molecule is already unbound
			return


		elif molecule.attr('_chromBoundToFork'):
			# Molecule is bound to a fork
			(strandIndex, position, directionBool, extentForward, extentReverse
				) = molecule.attrs('_chromStrand', '_chromPosition',
				'_chromDirection', '_chromExtentForward',
				'_chromExtentReverse')

			(childStrandA, childStrandB) = self._strandChildrenIndexes[strandIndex]

			(regionParent, regionChildA, regionChildB) = self._forkedRegions(
				position, directionBool, extentForward, extentReverse)

			self._array[strandIndex, regionParent] = self._empty
			self._array[childStrandA, regionChildA] = self._empty
			self._array[childStrandB, regionChildB] = self._empty

			molecule.attrIs(
				_chromBound = False,
				_chromBoundToFork = False
				)


		else:
			# Molecule is bound to a normal location (not a fork)
			(strandIndex, position, directionBool, extentForward, extentReverse
				) = molecule.attrs('_chromStrand', '_chromPosition',
				'_chromDirection', '_chromExtentForward',
				'_chromExtentReverse')

			region = self._region(position, directionBool, extentForward, extentReverse)

			self._array[strandIndex, region] = self._empty

			molecule.attrIs(_chromBound = False)

	# TODO: insure these methods don't return _fork objects
	def moleculesBound(self):
		return self._objectsContainer.objects(_chromBound = ('==', True))


	def moleculesUnbound(self):
		return self._objectsContainer.objects(_chromBound = ('==', False))


	def moleculesBoundWithName(self, moleculeName):
		return self._objectsContainer.objectsWithName(moleculeName, 
			_chromBound = ('==', True))


	def moleculeBoundAtPosition(self, strand, position): # TODO: test
		strandIndex = self._strandNameToIndex[strand]
		
		index = self._array[strandIndex, position] - self._offset

		if index < 0:
			return None

		else:
			return self._objectsContainer._objectByGlobalIndex(index)


	def moleculeBoundOnFork(self, fork): # TODO: test
		forkStrand, forkPosition, forkDirection = fork.attrs('_chromStrand',
			'_chromPosition', '_chromDirection')

		if forkDirection: # == (-)
			position = forkPosition - 1

		else: # == (+)
			position = forkPosition + 1
		
		index = self._array[forkStrand, position] - self._offset

		if index < 0: # anything < 0 must refer to a special value
			return None

		else:
			return self._objectsContainer._objectByGlobalIndex(index)


	def moleculesBoundOverExtent(self, strand, position, direction, extentForward, extentReverse):
		strandIndex = self._strandNameToIndex[strand]
		directionBool = self._directionCharToBool[direction]

		region = self._region(position, directionBool, extentForward, extentReverse)

		indexes = np.setdiff1d(self._array[strandIndex, region],
			self._specialValues) - self._offset

		return self._objectsContainer._objectsByGlobalIndex(indexes)


	def moleculesBoundNearMolecule(self, molecule, extentForward, extentReverse):
		# TODO: decide how to handle this for molecules on/near forks (tempted to just ignore or raise)
		raise NotImplementedError()


	def moleculesBoundNearFork(self, fork, extentForward, extentReverse):
		raise NotImplementedError()


	def divideRegion(self, strand, start, stop):
		# Break an empty, continuous region of a strand into two child strands,
		# and return the two "fork" objects

		# NOTE: start to stop is inclusive, unlike "range"!
		strandParent = self._strandNameToIndex[strand]
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
			_chromStrand = strandParent,
			_chromPosition = start,
			_chromDirection = 1, # NOTE: I've chosen the convention that the fork "direction" is towards the parent (i.e. the natural direction of extension)
			)

		forkStop = self._objectsContainer.objectNew(
			'_fork',
			_chromStrand = strandParent,
			_chromPosition = stop,
			_chromDirection = 0,
			)

		self._array[strandParent, start] = forkStart.attr('_globalIndex') + self._offset
		self._array[strandParent, stop] = forkStop.attr('_globalIndex') + self._offset

		return forkStart, forkStop


	def forks(self):
		# Return a set of forks
		return self._objectsContainer.objectsWithName('_fork')


	def forkExtend(self, fork, extent):
		# Move a fork along its parent strand (shortening the parent and 
		# elongating the children), and return the new position

		forkStrand, forkPosition, forkDirection = fork.attrs('_chromStrand',
			'_chromPosition', '_chromDirection')

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
		fork.attrIs(_chromPosition = newPosition)
		self._array[forkStrand, newPosition] = fork.attr('_globalIndex') + self._offset

		return newPosition


	def forksCombine(self, fork1, fork2):
		# Combine two forks on the same strand, splitting the chromosome
		raise NotImplementedError()


	def rootStrand(self):
		# Return the character used to identify the root strand
		return self._rootChar


	def regionsNearForks(self, extentForward, extentReverse, includeMoleculesOnEnds):
		forkedRegions = []

		forks = self.forks()

		forbiddenValues = [fork.attr('_globalIndex') + self._offset for fork in forks] + [self._inactive]

		for fork in self.forks():
			forkStrand, forkPosition, forkDirection = fork.attrs(
				'_chromStrand', '_chromPosition', '_chromDirection')

			regionParent, regionChildA, regionChildB = self._forkedRegions(
				forkPosition, forkDirection, extentForward, extentReverse)

			childStrandA, childStrandB = self._strandChildrenIndexes[forkStrand]

			# Trim any inactive/forked portions of the regions
			## on parent
			valuesParent, indexesParent = np.unique(self._array[forkStrand, regionParent], return_index = True)

			try:
				trimParentTo = np.min([
					index for value, index in zip(valuesParent, indexesParent) if value in forbiddenValues
					])

			except ValueError:
				pass

			else:
				regionParent = regionParent[:trimParentTo]


			## on child A
			valuesChildA, indexesChildA = np.unique(self._array[childStrandB, regionChildA], return_index = True)

			try:
				trimChildATo = np.min([
					index for value, index in zip(valuesChildA, indexesChildA) if value in forbiddenValues
					])

			except ValueError:
				pass

			else:
				regionChildA = regionChildA[:trimParentTo]


			## on child B
			valuesChildB, indexesChildB = np.unique(self._array[childStrandB, regionChildB], return_index = True)

			try:
				trimChildBTo = np.min([
					index for value, index in zip(valuesChildB, indexesChildB) if value in forbiddenValues
					])

			except ValueError:
				pass

			else:
				regionChildB = regionChildB[:trimParentTo]

			# Check to see whether the ends of the regions overlap with a 
			# molecule, and update accordingly

			forkMolecule = self.moleculeBoundOnFork(fork)

			## for parent
			parentEndMolecule = self.moleculeBoundAtPosition(
				self._strandNames[forkStrand], regionParent[-1])

			if parentEndMolecule is not None:
				footprint = self.moleculeFootprint(parentEndMolecule)

				if np.setdiff1d(footprint, regionParent).any():
					if forkMolecule is not None and parentEndMolecule == forkMolecule:
						raise ChrosomeContainerException('Something terrible happened')

					if not includeMoleculesOnEnds or parentEndMolecule.attr('_chromBoundToFork'):
						regionParent = np.array([i for i in regionParent if i not in footprint])

					else:
						regionParent = np.lib.arraysetops.union1d(footprint, regionParent)

						if forkDirection: # == (-)
							regionParent = np.flipud(regionParent)

			## for child A
			childAEndMolecule = self.moleculeBoundAtPosition(
				self._strandNames[childStrandA], regionChildA[-1])

			if childAEndMolecule is not None:
				footprint = self.moleculeFootprint(childAEndMolecule)

				if np.setdiff1d(footprint, regionChildA).any():
					if forkMolecule is not None and childAEndMolecule == forkMolecule:
						raise ChrosomeContainerException('Something terrible happened')

					if not includeMoleculesOnEnds or childAEndMolecule.attr('_chromBoundToFork'):
						regionChildA = np.array([i for i in regionChildA if i not in footprint])

					else:
						regionChildA = np.lib.arraysetops.union1d(footprint, regionChildA)

						if not forkDirection: # == (+)
							regionChildA = np.flipud(regionChildA)

			## for child B
			childBEndMolecule = self.moleculeBoundAtPosition(
				self._strandNames[childStrandB], regionChildB[-1])

			if childBEndMolecule is not None:
				footprint = self.moleculeFootprint(childBEndMolecule)

				if np.setdiff1d(footprint, regionChildB).any():
					if forkMolecule is not None and childBEndMolecule == forkMolecule:
						raise ChrosomeContainerException('Something terrible happened')

					if not includeMoleculesOnEnds or childBEndMolecule.attr('_chromBoundToFork'):
						regionChildB = np.array([i for i in regionChildB if i not in footprint])

					else:
						regionChildB = np.lib.arraysetops.union1d(footprint, regionChildB)

						if not forkDirection: # == (+)
							regionChildB = np.flipud(regionChildB)

			# Append region
			forkedRegions.append((forkStrand, regionParent[0],
				regionParent[-1], forkPosition, regionChildA[-1], regionChildB[-1]))

		# TODO: return a forked regions set object

		return _ForkedRegionSet(forkedRegions)


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

class _ForkedRegionSet(object):
	def __init__(self, forkedRegions):
		self._forkedRegions = np.array(forkedRegions,[
			('strand', np.int64),
			('parentStart', np.int64),
			('parentStop', np.int64),
			('childStart', np.int64),
			('childAStop', np.int64),
			('childBStop', np.int64),
			])


	def __str__(self):
		return str(self._forkedRegions)


		

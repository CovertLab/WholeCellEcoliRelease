
from __future__ import division

import numpy as np
import tables

import wholecell.states.state
import wholecell.utils.unique_objects_container

N_BASES = 5000000 # TODO: from kb
STRAND_MULTIPLICITY = 3
MOLECULE_WIDTH = 50 # TODO: from kb

MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		},
	'DNA polymerase':{
		},
	}

class SequenceBoundMolecules(object):
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
		# '_sequenceBoundToFork':'bool', # special property for molecules bound on a forked region
		}

	# Single-character names for denoting strand identity
	_rootChar = 'R' 
	_childLeadingChar = 'A'
	_childLaggingChar = 'B'

	# Conversion key for directions
	_positiveChar = '+'
	_negativeChar = '-'
	_directionCharToBool = {_positiveChar:False, _negativeChar:True}
	_directionBoolToChar = [_positiveChar, _negativeChar]
	
	def __init__(self):
		self._length = N_BASES

		self._strandMultiplicity = STRAND_MULTIPLICITY
		self._buildStrandConnectivity()
		
		self._array = np.zeros((self._nStrands, self._length), dtype = np.int32) # TODO: choose best dtype based on array size

		self._array[0, :] = self._empty # Root strand is always active

		moleculeAttributes = self._defaultObjectContainerObjects
		for moleculeName, attributes in MOLECULE_ATTRIBUTES.viewitems():
			moleculeAttributes[moleculeName] = attributes.copy()
			moleculeAttributes[moleculeName].update(self._defaultObjectContainerAttributes)

		self._objectsContainer =  wholecell.utils.unique_objects_container.UniqueObjectsContainer(
			moleculeAttributes)


	def _buildStrandConnectivity(self):
		self._strandNames = []
		self._strandNames.append(self._rootChar)

		parentNames = self._strandNames
		for i in xrange(self._strandMultiplicity - 1):
			childNames = []

			for parentName in parentNames:
				childNames.append(parentName + self._childLeadingChar)
				childNames.append(parentName + self._childLaggingChar)

			self._strandNames.extend(childNames)
			parentNames = childNames
		
		self._nStrands = len(self._strandNames)

		self._strandNameToIndex = {name:ind for ind, name in enumerate(self._strandNames)}
		self._strandChildrenIndexes = [
			(self._strandNameToIndex[name + self._childLeadingChar], self._strandNameToIndex[name + self._childLaggingChar])
			if self._strandNameToIndex.has_key(name + self._childLeadingChar) else None
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
		self._objectsContainer.objectDel(molecule)


	def moleculeLocationIs(self, molecule, strand, position, direction, extentForward, extentReverse):
		self.moleculeLocationIsUnbound(molecule)

		strandIndex = self._strandNameToIndex[strand]
		directionBool = self._directionCharToBool[direction]

		molecule.attrIs('_sequenceBound', True) # TODO: attrsAre method
		molecule.attrIs('_sequenceStrand', strandIndex)
		molecule.attrIs('_sequencePosition', position)
		molecule.attrIs('_sequenceDirection', directionBool)
		molecule.attrIs('_sequenceExtentForward', extentForward)
		molecule.attrIs('_sequenceExtentReverse', extentReverse)

		extentPositive, extentNegative = self._extentRelativeToAbsolute(
			extentForward, extentReverse, directionBool)

		region = self._array[strandIndex, (position - extentNegative):(position + extentPositive)]

		assert (region == self._empty).all(), 'Attempted to place a molecule in a non-empty region'

		region[:] = molecule.attr('_globalIndex') + self._offset


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

		else:
			strandIndex = molecule.attr('_sequenceStrand')
			position = molecule.attr('_sequencePosition')
			directionBool = molecule.attr('_sequenceDirection')
			extentForward = molecule.attr('_sequenceExtentForward')
			extentReverse = molecule.attr('_sequenceExtentReverse')

			extentPositive, extentNegative = self._extentRelativeToAbsolute(
				extentForward, extentReverse, directionBool)

			self._array[
				strandIndex,
				(position - extentNegative):(position + extentPositive)
				] = self._empty


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

				indexes = np.setdiff1d(
					self._array[strandIndex, (position - extentNegative):(position + extentPositive)],
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
		raise NotImplementedError()


	def divideRegion(self, strandName, start, stop):
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
			region = np.s_[forkPosition+1:forkPosition+extent+1]
			newPosition = forkPosition+extent

		else:
			region = np.s_[forkPosition-extent:forkPosition]
			newPosition = forkPosition-extent

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


	def moleculeLocationIsFork(self, molecule, fork):
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

# TODO: design views, queries, requests, etc

class Chromosome(wholecell.states.state.State):

	def __init__(self, *args, **kwargs):
		self.meta = {
			'id':'Chromosome',
			'name':'Chromosome',
			'dynamics':[],
			'units':{}
			}

		# self.time = None

		self.container = None

		super(Chromosome, self).__init__(*args, **kwargs)

	
	def initialize(self, sim, kb):
		super(Chromosome, self).initialize(sim, kb)

		self.container = SequenceBoundMolecules()


	def calcInitialConditions(self):
		# Add replication fork

		forkStartPos = 20000
		forkStopPos = 30000 - 1

		print (self.container._array[0, :] == self.container._empty).sum()

		(forkStart, forkStop) = self.container.divideRegion('R', forkStartPos, forkStopPos)

		print (self.container._array[0, :] == self.container._empty).sum(), forkStopPos - forkStartPos + 1
		print (self.container._array[1, :] == self.container._empty).sum()
		print (self.container._array[2, :] == self.container._empty).sum()

		self.container.forkExtend(forkStart, 20)

		print (self.container._array[0, :] == self.container._empty).sum()
		print (self.container._array[1, :] == self.container._empty).sum()
		print (self.container._array[2, :] == self.container._empty).sum()

		self.container.forkExtend(forkStop, 40)

		print (self.container._array[0, :] == self.container._empty).sum()
		print (self.container._array[1, :] == self.container._empty).sum()
		print (self.container._array[2, :] == self.container._empty).sum()
		
		# molecule = self.container.moleculeNew('DNA polymerase')

		# assert self.container.moleculeLocation(molecule) is None

		# self.container.moleculeLocationIs(
		# 	molecule,
		# 	'R', 100, '+',
		# 	50, 50
		# 	)

		# self.container.moleculeLocationIs(
		# 	molecule,
		# 	'R', 200, '+',
		# 	50, 50
		# 	)

		# assert self.container.moleculeLocation(molecule) == ('R', 100, '+', 50, 50)

		# assert self.container.moleculesBound() == {molecule}
		# assert self.container.moleculesBound('DNA polymerase') == {molecule}
		# assert self.container.moleculesBound('RNA polymerase') == set()

		# assert self.container.moleculesBound(
		# 	strand = 'R',
		# 	position = 100,
		# 	) == {molecule}

		# assert self.container.moleculesBound(
		# 	strand = 'R',
		# 	position = 140,
		# 	) == {molecule}

		# assert self.container.moleculesBound(
		# 	strand = 'R',
		# 	position = 160,
		# 	) == set()

		# assert self.container.moleculesBound(
		# 	moleculeName = 'DNA polymerase',
		# 	strand = 'R',
		# 	position = 100,
		# 	) == {molecule}

		# assert self.container.moleculesBound(
		# 	moleculeName = 'RNA polymerase',
		# 	strand = 'R',
		# 	position = 100,
		# 	) == set()

		# assert self.container.moleculesBound(
		# 	moleculeName = 'DNA polymerase',
		# 	strand = 'R',
		# 	position = 100,
		# 	direction = '+',
		# 	extentForward = 10,
		# 	extentReverse = 10,
		# 	) == {molecule}

		# assert self.container.moleculesBound(
		# 	moleculeName = 'DNA polymerase',
		# 	strand = 'R',
		# 	position = 100,
		# 	direction = '-',
		# 	extentForward = 10,
		# 	extentReverse = 10,
		# 	) == {molecule}

		# assert self.container.moleculesBound(
		# 	moleculeName = 'DNA polymerase',
		# 	strand = 'R',
		# 	position = 160,
		# 	direction = '-',
		# 	extentForward = 20,
		# 	extentReverse = 0,
		# 	) == {molecule}

		# assert self.container.moleculesBound(
		# 	moleculeName = 'DNA polymerase',
		# 	strand = 'R',
		# 	position = 160,
		# 	direction = '+',
		# 	extentForward = 20,
		# 	extentReverse = 0,
		# 	) == set()

		import ipdb; ipdb.set_trace()


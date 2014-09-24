"""
unqiue_objects_container.py

UniqueObjectsContainer is a object that tracks the attributes of unique
objects, which are typically molecules.  It supports saving and loading by 
appending entries in a structured array to tables with the same sets of fields.

The UniqueObjectsContainer uses _UniqueObject objects to present a clean 
interface to a specific molecule"s attributes.
"""

from __future__ import division

from copy import deepcopy
import warnings
from itertools import izip
from functools import partial

import numpy as np
import tables

# TODO: object transfer between UniqueObjectsContainer instances
# TODO: unique id for each object based on
#	hash(time, collectionIndex, arrayIndex, containerName) at creation time

_MAX_ID_SIZE = 40 # max length of the unique id assigned to objects

class UniqueObjectsContainerException(Exception):
	pass


class UniqueObjectsContainer(object):
	"""
	UniqueObjectsContainer

	Essentially a dict of structured arrays, where the structured array
	fields are attributes of the unique molecules.  Used for the unique
	molecules state and partitions.
	"""

	# State descriptions
	_entryInactive = 0 # a clear entry
	_entryActive = 1 # an entry that is in use

	_defaultSpecification = {
		"_entryState":np.int64, # see state descriptions above
		"_globalIndex":np.int64, # index in the _globalReference array (collection)
		"_timeStep":np.int64, # current time (important for saving) # TODO: handle, set in other classes
		"_uniqueId":"{}str".format(_MAX_ID_SIZE) # unique ID assigned to each object
		}

	_globalReferenceDtype = {
		"_entryState":np.int64, # see state descriptions above
		"_collectionIndex":np.int64,
		"_objectIndex":np.int64,
		"_timeStep":np.int64,
		}

	_fractionExtendEntries = 0.1 # fractional rate to increase number of entries in the structured array (collection)

	_queryOperations = {
		">":np.greater,
		">=":np.greater_equal,
		"<":np.less,
		"<=":np.less_equal,
		"==":np.equal,
		"!=":np.not_equal,
		"in":np.lib.arraysetops.in1d,
		"not in":partial(np.lib.arraysetops.in1d, invert = True)
		}

	def __init__(self, specifications):
		self._timeStep = 0

		self._names = [] # sorted list of object names

		self._collections = [] # ordered list of collections (which are arrays)
		self._nameToIndexMapping = {} # collectionName:index of associated structured array

		self._tableNames = [] # collection index:table name

		self._specifications = deepcopy(specifications) # collectionName:{attributeName:type}

		self._names = sorted(self._specifications.keys())

		defaultSpecKeys = self._defaultSpecification.viewkeys()

		# Add the attributes used internally
		for name, specification in self._specifications.viewitems():
			# Make sure there is no overlap in attribute names
			invalidAttributeNames = (specification.viewkeys() & defaultSpecKeys)
			if invalidAttributeNames:
				raise UniqueObjectsContainerException(
					"Invalid attribute names in specification for {}: {}".format(
						name,
						', '.join(invalidAttributeNames)
						)
					)

			specification.update(self._defaultSpecification)

		for collectionIndex, collectionName in enumerate(self._names):
			specification = self._specifications[collectionName]

			# Create the collection (structured array)
			newArray = np.zeros( # start out empty
				1,
				dtype = [
					(attrName, attrType)
					for attrName, attrType in specification.viewitems()
					]
				)

			# Create references to collections
			self._collections.append(newArray)
			self._nameToIndexMapping[collectionName] = collectionIndex

			# Give the tables accessible names
			self._tableNames.append(collectionName.replace(" ", "_"))

		# Create an array which handles global references to objects in all collections
		self._globalReference = np.zeros(
			1,
			dtype = [
				(attrName, attrType)
				for attrName, attrType in self._globalReferenceDtype.viewitems()
				]
			)


	def _getFreeIndexes(self, collectionIndex, nObjects):
		# Returns indexes of unoccupied entries, extending the arrays when neccesary

		collection = self._collections[collectionIndex]

		freeCollectionIndexes = np.where(
			collection["_entryState"] == self._entryInactive
			)[0]

		nFreeCollectionIndexes = freeCollectionIndexes.size

		if nFreeCollectionIndexes < nObjects:
			oldSize = collection.size
			nNewEntries = max(
				np.int64(oldSize * self._fractionExtendEntries),
				nObjects - nFreeCollectionIndexes
				)

			self._collections[collectionIndex] = np.append(
				collection,
				np.zeros(nNewEntries, dtype = collection.dtype)
				)

			freeCollectionIndexes = np.concatenate((
				freeCollectionIndexes,
				oldSize + np.arange(nNewEntries)
				))

		freeGlobalIndexes = np.where(
			self._globalReference["_entryState"] == self._entryInactive
			)[0]

		nFreeGlobalIndexes = freeGlobalIndexes.size

		if nFreeGlobalIndexes < nObjects:
			oldSize = self._globalReference.size
			nNewEntries = max(
				np.int64(oldSize * self._fractionExtendEntries),
				nObjects - nFreeGlobalIndexes
				)

			self._globalReference = np.append(
				self._globalReference,
				np.zeros(nNewEntries, dtype = self._globalReference.dtype)
				)

			freeGlobalIndexes = np.concatenate((
				freeGlobalIndexes,
				oldSize + np.arange(nNewEntries)
				))

		self._setCollectionTimes()

		return freeCollectionIndexes[:nObjects], freeGlobalIndexes[:nObjects]


	def objectsNew(self, collectionName, nObjects, **attributes):
		collectionIndex = self._nameToIndexMapping[collectionName]
		objectIndexes, globalIndexes = self._getFreeIndexes(collectionIndex, nObjects)

		collection = self._collections[collectionIndex]

		uniqueObjectIds = [
			"{}{}-{}".format(collectionName, self._timeStep, objectIndex)
			for objectIndex in objectIndexes
			]

		maxObjectIdLength = (
			len(collectionName)
			+ np.floor(np.log10(max(self._timeStep, 1))) + 1
			+ np.floor(np.log10(max(objectIndexes.max(), 1)))+ 1
			+ 1
			)

		if maxObjectIdLength > _MAX_ID_SIZE:
			warnings.warn("Maximum allowable ID size exceeded")

		collection["_entryState"][objectIndexes] = self._entryActive
		collection["_globalIndex"][objectIndexes] = globalIndexes
		collection["_uniqueId"][objectIndexes] = uniqueObjectIds

		for attrName, attrValue in attributes.viewitems():
			collection[attrName][objectIndexes] = attrValue

		self._globalReference["_entryState"][globalIndexes] = self._entryActive
		self._globalReference["_collectionIndex"][globalIndexes] = collectionIndex
		self._globalReference["_objectIndex"][globalIndexes] = objectIndexes

		return _UniqueObjectSet(self, globalIndexes)


	def objectNew(self, collectionName, **attributes):
		(molecule,) = self.objectsNew(collectionName, 1, **attributes) # NOTE: tuple unpacking

		return molecule


	def objectsDel(self, objects):
		for obj in objects:
			self.objectDel(obj)


	def objectDel(self, obj):
		# TODO: profile this in processes that delete lots of objects
		collection = self._collections[obj._collectionIndex]
		collection[obj._objectIndex].fill(0)

		self._globalReference[obj._globalIndex].fill(0)


	def objects(self, **operations):
		# Return all objects, optionally evaluating a query on !!every!! molecule (generally not what you want to do)
		if operations:
			results = []

			for collectionIndex in xrange(len(self._collections)):
				results.append(self._queryObjects(collectionIndex, **operations))

			return _UniqueObjectSet(self, np.concatenate([
				self._collections[collectionIndex]["_globalIndex"][result]
				for collectionIndex, result in enumerate(results)
				]))

		else:
			return _UniqueObjectSet(self,
				np.where(self._globalReference["_entryState"] == self._entryActive)[0]
				)


	def objectsInCollection(self, collectionName, **operations):
		# Return all objects belonging to a collection and that optionally satisfy a set of attribute queries
		collectionIndex = self._nameToIndexMapping[collectionName]

		result = self._queryObjects(collectionIndex, **operations)

		return _UniqueObjectSet(self,
			self._collections[collectionIndex]["_globalIndex"][result]
			)


	def objectsInCollections(self, collectionNames, **operations):
		# Return all objects belonging to a set of collections that optionally satisfy a set of attribute queries

		collectionIndexes = [self._nameToIndexMapping[collectionName] for collectionName in collectionNames]
		results = []

		for collectionIndex in collectionIndexes:
			results.append(self._queryObjects(collectionIndex, **operations))

		return _UniqueObjectSet(self, np.concatenate([
			self._collections[collectionIndex]["_globalIndex"][result]
			for collectionIndex, result in izip(collectionIndexes, results)
			]))


	def _queryObjects(self, collectionIndex, **operations):
		# Performs a series of comparison operations on a collection, and 
		# returns a boolean-value array that is True where every comparison was True
		operations["_entryState"] = ("==", self._entryActive)
		collection = self._collections[collectionIndex]

		return reduce(
			np.logical_and,
			(
				self._queryOperations[operator](
					collection[attrName],
					queryValue
					)
				for attrName, (operator, queryValue) in operations.viewitems()
			)
		)


	def objectsByGlobalIndex(self, globalIndexes):
		return _UniqueObjectSet(self, globalIndexes)


	def objectByGlobalIndex(self, globalIndex):
		return _UniqueObject(self, globalIndex)


	def objectByGlobalIndexDel(self, globalIndex):
		self.objectsByGlobalIndexDel(np.array(globalIndex))


	def timeStepIs(self, timeStep):
		self._timeStep = timeStep

		self._setCollectionTimes()


	def _setCollectionTimes(self):
		for collection in self._collections:
			collection["_timeStep"] = self._timeStep

		self._globalReference["_timeStep"] = self._timeStep


	def __eq__(self, other):
		return np.all(
			(selfCollection == otherCollection).all()
			for (selfCollection, otherCollection) in izip(self._collections, other._collections)
			) and np.all(self._globalReference == other._globalReference)


	def pytablesCreate(self, h5file):
		for collectionIndex, collection in enumerate(self._collections):
			h5file.create_table(
				h5file.root,
				self._tableNames[collectionIndex],
				collection.dtype,
				title = self._names[collectionIndex],
				filters = tables.Filters(complevel = 9, complib = "zlib")
				)

		h5file.create_table(
			h5file.root,
			"_globalReference",
			self._globalReference.dtype,
			title = "_globalReference",
			filters = tables.Filters(complevel = 9, complib = "zlib")
			)

		h5file.root._globalReference.attrs.collectionNames = self._names


	def pytablesAppend(self, h5file):
		for collectionIndex, collection in enumerate(self._collections):
			entryTable = h5file.get_node("/", self._tableNames[collectionIndex])

			entryTable.append(collection)

			entryTable.flush()

		globalTable = h5file.get_node("/", "_globalReference")

		globalTable.append(self._globalReference)

		globalTable.flush()


	def pytablesLoad(self, h5file, timePoint):
		for collectionIndex, tableName in enumerate(self._tableNames):
			entryTable = h5file.get_node("/", tableName)

			entries = entryTable[entryTable.col("_timeStep") == timePoint]

			self._collections[collectionIndex] = entries

		globalTable = h5file.get_node("/", "_globalReference")

		globalEntries = globalTable[globalTable.col("_timeStep") == timePoint]

		self._globalReference = globalEntries


class _UniqueObject(object):
	"""
	_UniqueObject

	A wrapper around a row in a container, refering to a specific unique object.
	Primarily used as a way to manipulate individual molecules with a python
	object-like interface.
	"""
	
	__slots__ = ("_container", "_globalIndex", "_collectionIndex", "_objectIndex")


	def __init__(self, container, globalIndex):
		self._container = container
		self._globalIndex = globalIndex
		self._collectionIndex = container._globalReference[globalIndex]["_collectionIndex"]
		self._objectIndex = container._globalReference[globalIndex]["_objectIndex"]


	def name(self):
		return self._container._collectionNames[self._collectionIndex]


	def uniqueId(self):
		return self.attr("_uniqueId")


	def attr(self, attribute):
		entry = self._container._collections[self._collectionIndex][self._objectIndex]
		
		if entry["_entryState"] == self._container._entryInactive:
			raise UniqueObjectsContainerException("Attempted to access an inactive object.")

		if isinstance(entry[attribute], np.ndarray):
			# Prevent making a view instead of a copy
			return entry[attribute].copy()

		else:
			return entry[attribute]


	def attrs(self, *attributes):
		entry = self._container._collections[self._collectionIndex][self._objectIndex]
		
		if entry["_entryState"] == self._container._entryInactive:
			raise UniqueObjectsContainerException("Attempted to access an inactive object.")
		
		# See note in .attr
		return tuple(
			entry[attribute].copy() if isinstance(entry[attribute], np.ndarray) else entry[attribute]
			for attribute in attributes
			)


	def attrIs(self, **attributes):
		entry = self._container._collections[self._collectionIndex][self._objectIndex]
		
		if entry["_entryState"] == self._container._entryInactive:
			raise UniqueObjectsContainerException("Attempted to access an inactive object.")

		for attribute, value in attributes.viewitems():
			if isinstance(entry[attribute], np.ndarray):
				# Fix for the circumstance that the attribute is an ndarray - 
				# without the [:] assignment, only the first value will be 
				# assigned (probably a NumPy bug)
				entry[attribute][:] = value

			else:
				entry[attribute] = value


	def __hash__(self):
		# TODO: replace with unique ID
		return hash((self._collectionIndex, self._objectIndex))


	def __eq__(self, other):
		if not isinstance(other, _UniqueObject):
			return False
			
		if not self._container is other._container:
			raise UniqueObjectsContainerException("Object comparisons across UniqueMoleculesContainer objects not supported.")

		return self._globalIndex == other._globalIndex


	def __ne__(self, other):
		return not self.__eq__(other)


class _UniqueObjectSet(object):
	"""
	_UniqueObjectSet

	A set of objects, stored internally by their global indexes.  Iterable and
	ordered.  Accessors allow for manipulating sets in lump.
	"""

	def __init__(self, container, globalIndexes):
		self._container = container
		self._globalIndexes = np.array(globalIndexes, np.int)


	def __contains__(self, uniqueObject):
		if not self._container is uniqueObject._container:
			raise UniqueObjectsContainerException("Object comparisons across UniqueMoleculesContainer objects not supported.")

		return uniqueObject.attr("_globalIndex") in self._globalIndexes


	def __iter__(self):
		return (_UniqueObject(self._container, globalIndex)
			for globalIndex in self._globalIndexes)


	def __len__(self):
		return self._globalIndexes.size


	def __eq__(self, other):
		if not self._container is other._container:
			raise UniqueObjectsContainerException("Object comparisons across UniqueMoleculesContainer objects not supported.")

		return (self._globalIndexes == other._globalIndexes).all()


	def __or__(self, other):
		if not self._container is other._container:
			raise UniqueObjectsContainerException("Object comparisons across UniqueMoleculesContainer objects not supported.")

		return _UniqueObjectSet(
			self._container,
			np.lib.arraysetops.union1d(self._globalIndexes, other._globalIndexes)
			)


	def __getitem__(self, index):
		return _UniqueObject(self._container, self._globalIndexes[index])


	def uniqueIds(self):
		return self.attr("_uniqueId")


	def attr(self, attribute):
		if self._globalIndexes.size == 0:
			raise UniqueObjectsContainerException("Object set is empty")

		if (self._container._globalReference["_entryState"][self._globalIndexes] == self._container._entryInactive).any():
			raise UniqueObjectsContainerException("One or more object was deleted from the set")

		# TODO: cache these properties? should be static
		collectionIndexes = self._container._globalReference["_collectionIndex"][self._globalIndexes]
		objectIndexes = self._container._globalReference["_objectIndex"][self._globalIndexes]

		uniqueColIndexes, inverse = np.unique(collectionIndexes, return_inverse = True)

		attributeDtype = self._container._collections[uniqueColIndexes[0]].dtype[attribute]

		values = np.zeros(
			self._globalIndexes.size,
			dtype = attributeDtype
			)

		for i, collectionIndex in enumerate(uniqueColIndexes):
			globalObjIndexes = np.where(inverse == i)
			objectIndexesInCollection = objectIndexes[globalObjIndexes]
			
			values[globalObjIndexes] = self._container._collections[collectionIndex][attribute][objectIndexesInCollection]

		return values


	def attrs(self, *attributes):
		return tuple(
			self.attr(attribute) for attribute in attributes
			)


	def attrsAsStructArray(self, *attributes):

		if self._globalIndexes.size == 0:
			raise UniqueObjectsContainerException("Object set is empty")

		if (self._container._globalReference["_entryState"][self._globalIndexes] == self._container._entryInactive).any():
			raise UniqueObjectsContainerException("One or more object was deleted from the set")

		# TODO: cache these properties? should be static
		collectionIndexes = self._container._globalReference["_collectionIndex"][self._globalIndexes]
		objectIndexes = self._container._globalReference["_objectIndex"][self._globalIndexes]

		uniqueColIndexes, inverse = np.unique(collectionIndexes, return_inverse = True)

		if len(attributes) == 0:
			attributes = self._container._collections[uniqueColIndexes[0]].dtype.names

		attributes = list(attributes)

		attributeDtypes = [
			(attribute, self._container._collections[uniqueColIndexes[0]].dtype[attribute])
			for attribute in attributes
			]

		values = np.zeros(
			self._globalIndexes.size,
			dtype = attributeDtypes
			)

		for i, collectionIndex in enumerate(uniqueColIndexes):
			globalObjIndexes = np.where(inverse == i)
			objectIndexesInCollection = objectIndexes[globalObjIndexes]
			
			values[globalObjIndexes] = self._container._collections[collectionIndex][attributes][objectIndexesInCollection]

		return values


	def attrIs(self, **attributes):
		if self._globalIndexes.size == 0:
			raise UniqueObjectsContainerException("Object set is empty")

		if (self._container._globalReference["_entryState"][self._globalIndexes] == self._container._entryInactive).any():
			raise UniqueObjectsContainerException("One or more object was deleted from the set")

		# TODO: cache these properties? should be static
		collectionIndexes = self._container._globalReference["_collectionIndex"][self._globalIndexes]
		objectIndexes = self._container._globalReference["_objectIndex"][self._globalIndexes]

		uniqueColIndexes, inverse = np.unique(collectionIndexes, return_inverse = True)

		for i, collectionIndex in enumerate(uniqueColIndexes):
			globalObjIndexes = np.where(inverse == i)
			objectIndexesInCollection = objectIndexes[globalObjIndexes]

			for attribute, values in attributes.viewitems():
				valuesAsArray = np.array(values, ndmin = 1)

				if valuesAsArray.shape[0] == 1: # is a singleton
					self._container._collections[collectionIndex][attribute][objectIndexesInCollection] = valuesAsArray

				else:
					self._container._collections[collectionIndex][attribute][objectIndexesInCollection] = valuesAsArray[globalObjIndexes]


	def delByIndexes(self, indexes):

		globalIndexes = self._globalIndexes[indexes]

		# TODO: cache these properties? should be static
		collectionIndexes = self._container._globalReference["_collectionIndex"][globalIndexes]
		objectIndexes = self._container._globalReference["_objectIndex"][globalIndexes]

		uniqueColIndexes, inverse = np.unique(collectionIndexes, return_inverse = True)

		for i, collectionIndex in enumerate(uniqueColIndexes):
			globalObjIndexes = np.where(inverse == i)
			objectIndexesInCollection = objectIndexes[globalObjIndexes]

			self._container._collections[collectionIndex][objectIndexesInCollection] = np.zeros(1, dtype=self._container._collections[collectionIndex].dtype)

		self._container._globalReference[globalIndexes] = np.zeros(1, dtype=self._container._globalReference.dtype)

	# TODO: set-like operations (union, intersection, etc.)


import wholecell.utils.linear_programming as lp

def _partition(objectRequestsArray, requestNumberVector, requestProcessArray, randomState):
	# Arguments:
	# objectRequestsArray: 2D bool array, (molecule)x(request)
	# requestNumberVector: number of molecules request, by request
	# requestProcessArray: 2D bool array, (request)x(process)
	# Returns:
	# partitionedMolecules: 2D bool array, (molecule)x(process)

	# TODO: full documentation/writeup, better docstring
	
	# Build matrix for optimization

	nObjects = objectRequestsArray.shape[0]
	nRequests = requestNumberVector.size
	nProcesses = requestProcessArray.shape[1]

	if nProcesses == 0:
		# Return nothing
		return np.zeros((nObjects, nProcesses), np.bool)

	# Make into structured array to condense the problem into unique rows
	objectRequestsStructured = objectRequestsArray.view(
		dtype = objectRequestsArray.dtype.descr * nRequests)

	uniqueEntriesStructured, mapping = np.unique(objectRequestsStructured,
		return_inverse = True)

	uniqueEntries = uniqueEntriesStructured.view((np.bool, nRequests))

	counts = np.bincount(mapping) # the number of each condensed molecule type

	nObjectTypes = counts.size

	# Some index mapping voodoo
	where0, where1 = np.where(uniqueEntries)

	nConnections = where0.size

	argsort = np.argsort(where1)

	moleculeToRequestConnections = np.zeros((nObjectTypes + nRequests,
		nConnections), np.int64)

	upperIndices = (where0, np.arange(where1.size)[argsort])
	lowerIndices = (nObjectTypes + where1[argsort], np.arange(where1.size))
	# End voodoo

	moleculeToRequestConnections[upperIndices] = -1
	moleculeToRequestConnections[lowerIndices] = 1

	# Create the matrix and fill in the values
	matrix = np.zeros(
		(nObjectTypes + nRequests + nProcesses,
			nObjectTypes + nConnections + 2*nProcesses),
		np.int64
		)

	# Molecule "boundary fluxes"
	matrix[:nObjectTypes, :nObjectTypes] = np.identity(nObjectTypes)

	# Flow from molecule type to request
	matrix[:nObjectTypes + nRequests,
		nObjectTypes:nObjectTypes+nConnections] = moleculeToRequestConnections

	# Flow from request to process
	matrix[nObjectTypes:nObjectTypes+nRequests,
		nObjectTypes+nConnections:nObjectTypes+nConnections+nProcesses][np.where(requestProcessArray)] = -requestNumberVector

	matrix[nObjectTypes + nRequests:,
		nObjectTypes+nConnections:nObjectTypes+nConnections+nProcesses] = np.identity(nProcesses)

	# Process "boundary fluxes"
	matrix[nObjectTypes + nRequests:,
		-nProcesses:] = -np.identity(nProcesses)

	# Create other linear programming parameters

	objective = np.zeros(matrix.shape[1], np.float)
	objective[-nProcesses:] = 1 # objective is to maximize process satisfaction
	# TODO: experiment with non-unity process weightings

	b = np.zeros(matrix.shape[0], np.float) # conservation law, i.e. b = 0 = Ax

	lowerBound = np.zeros(matrix.shape[1], np.float) # matrix is defined such that all values are >= 0

	upperBound = np.empty(matrix.shape[1], np.float)
	upperBound[:] = np.inf
	upperBound[:nObjectTypes] = counts # can use up to the total number of molecules
	upperBound[-nProcesses:] = 1 # processes can be up to 100% satisfied

	# Optimize

	solution = lp.linearProgramming(
		"maximize", objective,
		matrix.astype(np.float), b, # cvxopt requres floats
		lowerBound, upperBound,
		"S", "C", # no idea what these are supposed to do
		None # no options
		)[0].flatten()

	# Convert solution to amounts allocated to each process

	countsOfMoleculeTypeByConnection = -moleculeToRequestConnections[:nObjectTypes, :] * solution[nObjectTypes:nObjectTypes+nConnections]

	countsOfMoleculeTypeByConnection = np.floor(countsOfMoleculeTypeByConnection) # Round down to prevent oversampling

	countsOfMoleculeTypeByRequests = np.dot(
		countsOfMoleculeTypeByConnection,
		moleculeToRequestConnections[nObjectTypes:, :].T
		)

	countsOfMoleculeTypeByProcess = np.dot(
		countsOfMoleculeTypeByRequests,
		requestProcessArray
		)

	processOffsetsOfSelectedObjects = np.c_[
		np.zeros(nObjectTypes),
		np.cumsum(countsOfMoleculeTypeByProcess, axis = 1)
		].astype(np.int64)

	# TODO: find a way to eliminate the for-loops!
	partitionedMolecules = np.zeros((nObjects, nProcesses), np.bool)
	
	for moleculeIndex in np.arange(uniqueEntriesStructured.size):
		indexesOfSelectedObjects = np.where(moleculeIndex == mapping)[0]
		randomState.shuffle(indexesOfSelectedObjects)

		for processIndex in np.arange(nProcesses):
			start = processOffsetsOfSelectedObjects[moleculeIndex, processIndex]
			stop = processOffsetsOfSelectedObjects[moleculeIndex, processIndex + 1]
			selectedIndexes = indexesOfSelectedObjects[start : stop]

			partitionedMolecules[
				selectedIndexes,
				processIndex] = True
	
	return partitionedMolecules

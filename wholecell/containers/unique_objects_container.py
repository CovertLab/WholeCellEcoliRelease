'''
unqiue_objects_container.py

UniqueObjectsContainer is a object that tracks the attributes of unique
objects, which are typically molecules.  It supports saving and loading by 
appending entries in a structured array to tables with the same sets of fields.

The UniqueObjectsContainer uses _UniqueObject objects to present a clean 
interface to a specific molecule's attributes.
'''

from __future__ import division

import numpy as np
import tables

import wholecell.utils.linear_programming as lp

# State descriptions
_ENTRY_INACTIVE = 0 # a clear entry
_ENTRY_ACTIVE = 1 # an entry that is in use
_ENTRY_DELETED = 2 # an entry that was deleted and is waiting to be cleaned up

# TODO: reevaluate the assignment of private methods/attributes in this and similar classes

class UniqueObjectsContainerException(Exception):
	pass


class UniqueObjectsContainer(object):
	'''
	UniqueObjectsContainer

	Essentially a dict of structured arrays, where the structured array
	fields are attributes of the unique molecules.  Used for the unique
	molecules state and partitions.
	'''

	_defaultContainerAttributes = {
		'_entryState':'int64', # see state descriptions above
		'_globalIndex':'int64', # index in the _globalReference array (collection)
		'_time':'int64', # current time (important for saving)
		'_partitionedProcess':'uint32' # TODO: assign to every as default instead of in this fake KB
		}

	_defaultCollectionsSpec = {
		'_globalReference':{ # a table which contains reference to all molecules
			'_collectionIndex':'int64',
			'_objectIndex':'int64'
			}
		}

	_fractionExtendEntries = 0.1 # fractional rate to increase number of entries in the structured array (collection)

	_queryOperations = {
		'>':np.greater,
		'>=':np.greater_equal,
		'<':np.less,
		'<=':np.less_equal,
		'==':np.equal,
		'!=':np.not_equal
		}

	def __init__(self, collectionsSpec):
		self._time = 0

		self._collectionsSpec = {} # collectionName:{attributeName:type}

		self._collectionNames = [] # sorted list of object names

		self._collections = [] # ordered list of collections (which are arrays)
		self._collectionNameToIndexMapping = {} # collectionName:index of associated structured array

		self._tableNames = [] # collection index:table name

		self._collectionsSpec.update(collectionsSpec)
		self._collectionsSpec.update(self._defaultCollectionsSpec)

		self._collectionNames = sorted(self._collectionsSpec.keys())
		self._globalRefIndex = self._collectionNames.index('_globalReference')

		for collectionName, attributes in self._collectionsSpec.viewitems():
			# Add the attributes used internally
			attributes.update(self._defaultContainerAttributes)

		# Global references don't use global indexes
		del self._collectionsSpec['_globalReference']['_globalIndex']

		for collectionIndex, collectionName in enumerate(self._collectionNames):
			attributes = self._collectionsSpec[collectionName]

			# Create the collection (structured array)
			newArray = np.zeros(
				0, # start out empty
				dtype = [
					(attrName, attrType)
					for attrName, attrType in attributes.viewitems()
					]
				)

			# Create references to collections
			self._collections.append(newArray)
			self._collectionNameToIndexMapping[collectionName] = collectionIndex

			# Give the tables accessible names
			self._tableNames.append(collectionName.replace(' ', '_'))


	def objectsNew(self, collectionName, nMolecules, **attributes):
		attributes['_time'] = self._time

		# Create multiple objects of the same type and attribute values
		collectionIndex = self._collectionNameToIndexMapping[collectionName]
		objectIndexes = self._getFreeIndexes(collectionIndex, nMolecules)

		collection = self._collections[collectionIndex]

		collection['_entryState'][objectIndexes] = _ENTRY_ACTIVE

		for attrName, attrValue in attributes.viewitems():
			# NOTE: there is probably a non-loop solution to this, but the 'obvious' solution creates a copy instead of a view
			collection[attrName][objectIndexes] = attrValue

		globalIndexes = self._getFreeIndexes(self._globalRefIndex, nMolecules)

		# In global array, create reference pointing to object
		globalArray = self._collections[self._globalRefIndex]
		globalArray['_entryState'][globalIndexes] = _ENTRY_ACTIVE
		globalArray['_collectionIndex'][globalIndexes] = collectionIndex
		globalArray['_objectIndex'][globalIndexes] = objectIndexes
		globalArray['_time'][globalIndexes] = self._time

		# In collection, for each object, point to global reference
		collection['_globalIndex'][objectIndexes] = globalIndexes

		return _UniqueObjectSet(self, globalIndexes)


	def objectNew(self, collectionName, **attributes):
		(molecule,) = self.objectsNew(collectionName, 1, **attributes) # NOTE: tuple unpacking

		return molecule


	def _getFreeIndexes(self, collectionIndex, nMolecules):
		freeIndexes = np.where(
			self._collections[collectionIndex]['_entryState'] == _ENTRY_INACTIVE
			)[0]

		if freeIndexes.size < nMolecules:
			oldArray = self._collections[collectionIndex]
			oldSize = oldArray.size

			newSize = oldSize + max(int(oldSize * self._fractionExtendEntries), nMolecules)

			self._collections[collectionIndex] = np.zeros(
				newSize,
				dtype = oldArray.dtype
				)
			
			self._collections[collectionIndex][:oldSize] = oldArray

			freeIndexes = np.concatenate((freeIndexes, np.arange(oldSize, newSize)))

		return freeIndexes[:nMolecules]


	def objectsDel(self, objects):
		for obj in objects:
			self.objectDel(obj)


	def objectDel(self, obj):
		globalIndex = obj.attr('_globalIndex')

		self._collections[obj._collectionIndex][obj._objectIndex]['_entryState'] = _ENTRY_DELETED
		self._collections[self._globalRefIndex][globalIndex]['_entryState'] = _ENTRY_DELETED
		# TODO: Assign unique IDs and run _clearEntries() here


	def _clearEntries(self, collectionIndex, objectIndexes):
		collection = self._collections[collectionIndex]

		collection[objectIndexes] = np.zeros(
			1,
			dtype = collection.dtype
			)


	def objects(self, **operations):
		# Return all objects, optionally evaluating a query on !!every!! molecule (generally not what you want to do)
		if operations:
			collectionIndexes = set(xrange(len(self._collections)))
			collectionIndexes.remove(self._globalRefIndex)

			results = []

			for collectionIndex in collectionIndexes:
				results.append(self._queryObjects(collectionIndex,
					raiseOnMissingAttribute = False, **operations))

			return _UniqueObjectSet(self, np.r_[tuple(
				self._collections[collectionIndex]['_globalIndex'][result]
				for collectionIndex, result in zip(collectionIndexes, results)
				)])

		else:
			return _UniqueObjectSet(self,
				np.where(self._collections[self._globalRefIndex]['_entryState'] == _ENTRY_ACTIVE)[0]
				)


	def objectsInCollection(self, collectionName, **operations):
		# Return all objects belonging to a collection and that optionally satisfy a set of attribute queries
		collectionIndex = self._collectionNameToIndexMapping[collectionName]

		result = self._queryObjects(collectionIndex, **operations)

		return _UniqueObjectSet(self,
			self._collections[collectionIndex]['_globalIndex'][result]
			)


	def objectsInCollections(self, collectionNames, **operations):
		# Return all objects belonging to a set of collections that optionally satisfy a set of attribute queries

		collectionIndexes = [self._collectionNameToIndexMapping[collectionName] for collectionName in collectionNames]
		results = []

		for collectionIndex in collectionIndexes:
			results.append(self._queryObjects(collectionIndex, **operations))

		return _UniqueObjectSet(self, np.r_[tuple(
			self._collections[collectionIndex]['_globalIndex'][result]
			for collectionIndex, result in zip(collectionIndexes, results)
			)])


	def _queryObjects(self, collectionIndex, raiseOnMissingAttribute = True, **operations):
		operations['_entryState'] = ('==', _ENTRY_ACTIVE)
		collection = self._collections[collectionIndex]

		try:
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

		except ValueError:
			if not raiseOnMissingAttribute and (operations.viewkeys() - set(collection.dtype.names)):
				return np.zeros(collection.size, np.bool)

			else:
				raise


	def objectsByGlobalIndex(self, globalIndexes): # ASK JM: Make public?
		return _UniqueObjectSet(self, globalIndexes)


	def objectByGlobalIndex(self, globalIndex): # ASK JM: Make public?
		return _UniqueObject(self, globalIndex)


	def timeIs(self, time):
		self._time = time

		for collection in self._collections:
			collection['_time'] = time


	def flushDeleted(self):
		for collectionIndex, collection in enumerate(self._collections):
			self._clearEntries(
				collectionIndex,
				np.where(collection['_entryState'] == _ENTRY_DELETED)
				)


	def __eq__(self, other):
		return np.all(
			(selfCollection == otherCollection).all()
			for (selfCollection, otherCollection) in zip(self._collections, other._collections)
			)


	def pytablesCreate(self, h5file):
		for collectionIndex, collection in enumerate(self._collections):
			h5file.create_table(
				h5file.root,
				self._tableNames[collectionIndex],
				collection.dtype,
				title = self._collectionNames[collectionIndex],
				filters = tables.Filters(complevel = 9, complib = 'zlib')
				)


	def pytablesAppend(self, h5file):
		for collectionIndex, collection in enumerate(self._collections):
			entryTable = h5file.get_node('/', self._tableNames[collectionIndex])

			entryTable.append(collection)

			entryTable.flush()


	def pytablesLoad(self, h5file, timePoint):
		for collectionIndex, tableName in enumerate(self._tableNames):
			entryTable = h5file.get_node('/', tableName)

			entries = entryTable[entryTable.col('_time') == timePoint]

			self._collections[collectionIndex] = entries


	# TODO: compute mass
	# TODO: implement molecule transfer between containers


class _UniqueObject(object):
	'''
	_UniqueObject

	A wrapper around a row in a container, refering to a specific unique object.
	Primarily used as a way to manipulate individual molecules with a python
	object-like interface.
	'''
	
	__slots__ = ('_container', '_collectionIndex', '_objectIndex')


	def __init__(self, container, globalIndex):
		self._container = container
		self._collectionIndex = container._collections[container._globalRefIndex][globalIndex]['_collectionIndex']
		self._objectIndex = container._collections[container._globalRefIndex][globalIndex]['_objectIndex']


	def name(self):
		return self._container._collectionNames[self._collectionIndex]


	def attr(self, attribute):
		entry = self._container._collections[self._collectionIndex][self._objectIndex]
		
		if not entry['_entryState'] == _ENTRY_ACTIVE:
			raise UniqueObjectsContainerException('Attempted to access an inactive object.')

		if isinstance(entry[attribute], np.ndarray):
			# Prevent making a view instead of a copy
			return entry[attribute].copy()

		else:
			return entry[attribute]


	def attrs(self, *attributes):
		entry = self._container._collections[self._collectionIndex][self._objectIndex]
		
		if not entry['_entryState'] == _ENTRY_ACTIVE:
			raise UniqueObjectsContainerException('Attempted to access an inactive object.')
		
		# See note in .attr
		return tuple(
			entry[attribute].copy() if isinstance(entry[attribute], np.ndarray) else entry[attribute]
			for attribute in attributes
			)


	def attrIs(self, **attributes):
		entry = self._container._collections[self._collectionIndex][self._objectIndex]
		
		if not entry['_entryState'] == _ENTRY_ACTIVE:
			raise UniqueObjectsContainerException('Attempted to access an inactive object.')

		for attribute, value in attributes.viewitems():
			if isinstance(entry[attribute], np.ndarray):
				# Fix for the circumstance that the attribute is an ndarray - 
				# without the [:] assignment, only the first value will be 
				# assigned (probably a NumPy bug)
				entry[attribute][:] = value

			else:
				entry[attribute] = value


	def __hash__(self):
		return hash((self._collectionIndex, self._objectIndex))


	def __eq__(self, other):
		if not isinstance(other, _UniqueObject):
			return False
			
		if not self._container is other._container:
			raise UniqueObjectsContainerException('Object comparisons across UniqueMoleculesContainer objects not supported.')

		return self._collectionIndex == other._collectionIndex and self._objectIndex == other._objectIndex


	def __ne__(self, other):
		return not self.__eq__(other)


class _UniqueObjectSet(object):
	'''
	_UniqueObjectSet

	A set of objects, stored internally by their global indexes.  Iterable and
	ordered.  Accessors allow for manipulating sets in lump.
	'''

	# TODO: look into subclassing from collections.ViewKeys
	def __init__(self, container, globalIndexes):
		self._container = container
		self._globalIndexes = np.array(globalIndexes)


	def __contains__(self, uniqueObject):
		assert uniqueObject._container is self._container
		return uniqueObject.attr('_globalIndex') in self._globalIndexes


	def __iter__(self):
		return (_UniqueObject(self._container, globalIndex)
			for globalIndex in self._globalIndexes)


	def __len__(self):
		return self._globalIndexes.size


	def __eq__(self, other):
		if not self._container is other._container:
			raise UniqueObjectsContainerException('Object comparisons across UniqueMoleculesContainer objects not supported.')

		return (self._globalIndexes == other._globalIndexes).all()


	def __or__(self, other):
		assert self._container is other._container

		return _UniqueObjectSet(
			self._container,
			np.lib.arraysetops.union1d(self._globalIndexes, other._globalIndexes)
			)


	def attr(self, attribute):
		if self._globalIndexes.size == 0:
			return np.zeros(0) # NOTE: this does not enforce dtype; that is difficult!

		# TODO: cache these properties? should be static
		globalRef = self._container._collections[self._container._globalRefIndex]

		collectionIndexes = globalRef['_collectionIndex'][self._globalIndexes]
		objectIndexes = globalRef['_objectIndex'][self._globalIndexes]

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


	def attrIs(self, **attributes):
		if self._globalIndexes.size == 0:
			return

		# TODO: cache these properties? should be static
		globalRef = self._container._collections[self._container._globalRefIndex]

		collectionIndexes = globalRef['_collectionIndex'][self._globalIndexes]
		objectIndexes = globalRef['_objectIndex'][self._globalIndexes]

		uniqueColIndexes, inverse = np.unique(collectionIndexes, return_inverse = True)

		for i, collectionIndex in enumerate(uniqueColIndexes):
			globalObjIndexes = np.where(inverse == i)
			objectIndexesInCollection = objectIndexes[globalObjIndexes]

			for attribute, values in attributes.viewitems():
				valuesAsArray = np.array(values)

				if valuesAsArray.size == 1: # is a singleton
					self._container._collections[collectionIndex][attribute][objectIndexesInCollection] = valuesAsArray

				else:
					self._container._collections[collectionIndex][attribute][objectIndexesInCollection] = valuesAsArray[globalObjIndexes]

	# TODO: set-like operations (union, intersection, etc.)


def _partition(objectRequestsArray, requestNumberVector, requestProcessArray, randStream):
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
		'maximize', objective,
		matrix.astype(np.float), b, # cvxopt requres floats
		lowerBound, upperBound,
		'S', 'C', # no idea what these are supposed to do
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
		randStream.numpyShuffle(indexesOfSelectedObjects)

		for processIndex in np.arange(nProcesses):
			start = processOffsetsOfSelectedObjects[moleculeIndex, processIndex]
			stop = processOffsetsOfSelectedObjects[moleculeIndex, processIndex + 1]
			selectedIndexes = indexesOfSelectedObjects[start : stop]

			partitionedMolecules[
				selectedIndexes,
				processIndex] = True
	
	return partitionedMolecules

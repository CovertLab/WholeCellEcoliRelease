'''
unqiue_objects_container.py

UniqueObjectsContainer is a object that tracks the attributes of unique
objects, which are typically molecules.  It supports saving and loading by 
appending entries in a structured array to tables with the same sets of fields.

The UniqueObjectsContainer uses _UniqueObject objects to present a clean 
interface to a specific molecule's attributes.

The UniqueObjectsContainer also uses _Query objects to store and periodically
update queries which return sets of _UniqueObject that refer to molecules 
satisfying the query.

'''

from __future__ import division

import numpy as np
import tables

import wholecell.utils.linear_programming as lp

ENTRY_INACTIVE = 0 # a clear entry
ENTRY_ACTIVE = 1 # an entry that is in use
ENTRY_DELETED = 2 # an entry that was deleted and is waiting to be cleaned up

class UniqueObjectsContainer(object):
	'''
	UniqueObjectsContainer

	Essentially a dict of structured arrays, where the fields are attributes
	of the unique molecules.  Used for the unique molecules state and
	partitions.
	'''

	_defaultContainerAttributes = {
		'_entryState':'uint32', # see state descriptions above
		'_globalIndex':'uint32', # index in the _globalReference array
		'_time':'uint64', # current time (important for saving)
		# '_massDifference':'float64' # dynamic mass difference
		}

	_defaultObjects = {
		'_globalReference':{ # a table which contains reference to all molecules
			'_arrayIndex':'uint32',
			'_objectIndex':'uint32'
			}
		}

	_fractionExtendEntries = 0.1 # fractional rate to increase number of entries in the structured array

	_queryOperations = {
		'>':np.greater,
		'>=':np.greater_equal,
		'<':np.less,
		'<=':np.less_equal,
		'==':np.equal,
		'!=':np.not_equal
		}

	def __init__(self, objectAttributes):
		self._objectAttributes = {} # objectName:{attributeName:type}

		self._objectNames = [] # sorted list of object names

		self._arrays = [] # ordered list of arrays
		self._nameToArrayIndex = {} # objectName:index of associated structured array

		self._queries = [] # list of _Query objects
		self._tableNames = {} # objectName:table name

		self._objectAttributes.update(objectAttributes)
		self._objectAttributes.update(self._defaultObjects)

		self._objectNames = sorted(self._objectAttributes.keys())
		self._globalRefIndex = self._objectNames.index('_globalReference')

		for objectName, attributes in self._objectAttributes.viewitems():
			# Add the attributes used internally
			attributes.update(self._defaultContainerAttributes)

		# Global references don't use global indexes
		del self._objectAttributes['_globalReference']['_globalIndex']

		for arrayIndex, objectName in enumerate(self._objectNames):
			attributes = self._objectAttributes[objectName]

			# Create the structured array
			newArray = np.zeros(
				0, # start out empty
				dtype = [
					(attrName, attrType)
					for attrName, attrType in attributes.viewitems()
					]
				)

			# Create references to arrays
			self._arrays.append(newArray)
			self._nameToArrayIndex[objectName] = arrayIndex

			# Give the tables accessible names
			self._tableNames[objectName] = objectName.replace(' ', '_')

		# TODO: alternate constructor for copying to partitions


	def objectsNew(self, objectName, nMolecules, **attributes):
		arrayIndex = self._nameToArrayIndex[objectName]
		objectIndexes = self._getFreeIndexes(arrayIndex, nMolecules)

		array = self._arrays[arrayIndex]

		array['_entryState'][objectIndexes] = ENTRY_ACTIVE

		for attrName, attrValue in attributes.viewitems():
			# NOTE: there is probably a non-loop solution to this, but the 'obvious' solution creates a copy instead of a view
			array[attrName][objectIndexes] = attrValue

		globalIndexes = self._getFreeIndexes(self._globalRefIndex, nMolecules)
		globalArray = self._arrays[self._globalRefIndex]
		globalArray['_entryState'][globalIndexes] = ENTRY_ACTIVE
		globalArray['_arrayIndex'][globalIndexes] = arrayIndex
		globalArray['_objectIndex'][globalIndexes] = objectIndexes

		array['_globalIndex'][objectIndexes] = globalIndexes

		return self._objects(arrayIndex, objectIndexes)


	def objectNew(self, objectName, **attributes):
		(molecule,) = self.objectsNew(objectName, 1, **attributes) # NOTE: tuple unpacking

		return molecule


	def _getFreeIndexes(self, arrayIndex, nMolecules):
		freeIndexes = np.where(
			self._arrays[arrayIndex]['_entryState'] == ENTRY_INACTIVE
			)[0]

		if freeIndexes.size < nMolecules:
			oldEntries = self._arrays[arrayIndex]
			oldSize = oldEntries.size

			newSize = oldSize + max(int(oldSize * self._fractionExtendEntries), nMolecules)

			self._arrays[arrayIndex] = np.zeros(
				newSize,
				dtype = oldEntries.dtype
				)
			
			self._arrays[arrayIndex][:oldSize] = oldEntries

			freeIndexes = np.concatenate((freeIndexes, np.arange(oldSize, newSize)))


		return freeIndexes[:nMolecules]


	def objectsDel(self, objects):
		for obj in objects:
			self.objectDel(obj)


	def objectDel(self, obj):
		globalIndex = obj.attr('_globalIndex')

		self._arrays[obj._arrayIndex][obj._objectIndex]['_entryState'] = ENTRY_DELETED
		self._arrays[self._globalRefIndex][globalIndex]['_entryState'] = ENTRY_DELETED


	def _clearEntries(self, arrayIndex, objectIndexes):
		array = self._arrays[arrayIndex]

		array[objectIndexes] = np.zeros(
			1,
			dtype = array.dtype
			)


	def evaluateQuery(self, objectName, **operations): # TODO: allow for queries over all or a subset of molecules
		arrayIndex = self._nameToArrayIndex[objectName]

		return self._objects(
			arrayIndex,
			np.where(self._queryObjects(arrayIndex, **operations))[0]
			)
	

	def _queryObjects(self, arrayIndex, **operations):
		operations['_entryState'] = ('==', ENTRY_ACTIVE)
		array = self._arrays[arrayIndex]

		return reduce(
			np.logical_and,
			(
				self._queryOperations[operator](
					array[attrName],
					queryValue
					)
				for attrName, (operator, queryValue) in operations.viewitems()
			)
		)


	def updateQueries(self):
		for query in self._queries:
			query._objectIndexes = np.where(
				self._queryObjects(query._arrayIndex, **query._operations)
				)[0]


	def queryNew(self, objectName, **operations):
		arrayIndex = self._nameToArrayIndex[objectName]

		query = _Query(self, arrayIndex, **operations)

		self._queries.append(query)

		# TODO: return an old query object if the query is the same

		return query


	def objects(self, objectName):
		return self._objects(self._nameToArrayIndex[objectName])


	def _objects(self, arrayIndex, objectIndexes = None):
		return set(self._iterObjects(arrayIndex, objectIndexes)) # TODO: return a set-like object that creates the _UniqueObject instances as needed


	def _objectsByGlobalIndex(self, globalIndexes): # TODO: make global index ref the internal standard behavior?
		globalArray = self._arrays[self._globalRefIndex]

		return set(
			_UniqueObject(self, globalArray['_arrayIndex'][index], globalArray['_objectIndex'][index])
			for index in globalIndexes
			)


	def iterObjects(self, objectName):
		return self._iterObjects(self._nameToArrayIndex[objectName])


	def _iterObjects(self, arrayIndex, objectIndexes = None):
		if objectIndexes is None:
			objectIndexes = np.where(self._arrays[arrayIndex]['_entryState'] == ENTRY_ACTIVE)[0]

		return (_UniqueObject(self, arrayIndex, objectIndex) for objectIndex in objectIndexes)



	def _timeIs(self, time):
		for array in self._arrays:
			array['_time'] = time


	def _flushDeleted(self):
		for arrayIndex, array in enumerate(self._arrays):
			self._clearEntries(
				arrayIndex,
				np.where(array['_entryState'] == ENTRY_DELETED)
				)


	def __eq__(self, other):
		return all(
			(selfArray[selfArray['_entryState'] != ENTRY_INACTIVE] == otherArray[otherArray['_entryState'] != ENTRY_INACTIVE]).all()
			for (selfArray, otherArray) in zip(self._arrays, other._arrays)
			)


	def pytablesCreate(self, h5file):
		for arrayIndex, array in enumerate(self._arrays):
			h5file.create_table(
				h5file.root,
				self._tableNames[arrayIndex],
				array[self._savedAttributes[arrayIndex]].dtype,
				title = self._objectNames[arrayIndex],
				filters = tables.Filters(complevel = 9, complib = 'zlib')
				)

			h5file.create_table(
				h5file.root,
				self._tableNames[arrayIndex] + '_indexes',
				{'time':tables.UInt32Col(), 'index':tables.UInt32Col()},
				title = self._objectNames[arrayIndex] + ' indexes',
				filters = tables.Filters(complevel = 9, complib = 'zlib')
				)


	def pytablesAppend(self, h5file, time):
		for arrayIndex, array in enumerate(self._arrays):
			activeIndexes = np.where(array['_entryState'] != ENTRY_INACTIVE)[0]

			entryTable = h5file.get_node('/', self._tableNames[arrayIndex])

			entries = array[activeIndexes][self._savedAttributes[arrayIndex]]

			entryTable.append(entries)

			entryTable.flush()

			indexTable = h5file.get_node('/', self._tableNames[arrayIndex] + '_indexes')

			indexes = np.empty((activeIndexes.size, 2), [('time', np.uint32), ('indexes', np.uint32)])
			indexes['time'] = entries['time']
			indexes['indexes'] = activeIndexes

			indexTable.append(entries)

			indexTable.flush()


	def pytablesLoad(self, h5file, timePoint):
		for arrayIndex, tableName in enumerate(self._tableNames):
			entryTable = h5file.get_node('/', tableName)

			entries = entryTable[entryTable['_time'] == timePoint]

			indexTable = h5file.get_node('/', tableName + '_indexes')

			indexes = indexTable[indexTable['_time'] == timePoint]['indexes']

			self._arrays[arrayIndex] = np.array(
				indexes.max()+1,
				dtype = self._arrays[arrayIndex]
				)

			self._arrays[indexes] = entries


	# TODO: compute mass
	# TODO: implement molecule transfer between containers


class _Query(object):
	'''
	_Query

	A reference to a query that can be updated by the parent container and 
	inspected by a process.
	'''

	__slots__ = ('_container', '_arrayIndex', '_operations', '_objectIndexes')


	def __init__(self, container, arrayIndex, **operations):
		self._container = container
		self._arrayIndex = arrayIndex
		self._operations = operations
		self._objectIndexes = None


	def objects(self):
		return self._container._objects(self._arrayIndex, self._objectIndexes)


	def iterObjects(self):
		return self._container._iterObjects(self._arrayIndex, self._objectIndexes)


	# TODO: sampling functions?  i.e. get N molecules
	# TODO: subqueries?


class _UniqueObject(object):
	'''
	_UniqueObject

	A wrapper around a row in a container, refering to a specific unique object.
	Primarily used as a way to manipulate individual molecules with a python
	object-like interface.
	'''
	
	__slots__ = ('_container', '_arrayIndex', '_objectIndex')


	def __init__(self, container, arrayIndex, index):
		self._container = container
		self._arrayIndex = arrayIndex
		self._objectIndex = index


	def attr(self, attribute):
		entry = self._container._arrays[self._arrayIndex][self._objectIndex]
		
		if not entry['_entryState'] == ENTRY_ACTIVE:
			raise Exception('Attempted to access an inactive molecule.')

		return entry[attribute]


	def attrIs(self, attribute, value):
		entry = self._container._arrays[self._arrayIndex][self._objectIndex]
		
		if not entry['_entryState'] == ENTRY_ACTIVE:
			raise Exception('Attempted to access an inactive molecule.')

		entry[attribute] = value


	def __hash__(self):
		return hash((self._container, self._arrayIndex, self._objectIndex))


	def __eq__(self, other):
		assert self._container == other._container, 'Molecule comparisons across UniqueMoleculesContainer objects not supported.'
		return self._arrayIndex == other._arrayIndex and self._objectIndex == other._objectIndex

	# TODO: method to get name

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

	# Make into structured array to condense the problem into unique rows
	objectRequestsStructured = objectRequestsArray.view(
		dtype = [('', np.bool, nRequests)])

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
		nConnections), np.int)

	upperIndices = (where0, np.arange(where1.size)[argsort])
	lowerIndices = (nObjectTypes + where1[argsort], np.arange(where1.size))
	# End voodoo

	moleculeToRequestConnections[upperIndices] = -1
	moleculeToRequestConnections[lowerIndices] = 1

	# Create the matrix and fill in the values
	matrix = np.zeros(
		(nObjectTypes + nRequests + nProcesses,
			nObjectTypes + nConnections + 2*nProcesses),
		np.int
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

	unfixedCounts = -moleculeToRequestConnections[:nObjectTypes, :] * solution[nObjectTypes:nObjectTypes+nConnections]

	flooredCounts = np.floor(unfixedCounts) # Round down to prevent oversampling

	flooredProcessCounts = np.dot(
		np.dot(flooredCounts, moleculeToRequestConnections[nObjectTypes:, :].T),
		requestProcessArray
		)

	indexingRanges = np.c_[np.zeros(nObjectTypes), np.cumsum(flooredProcessCounts, 1)].astype(np.int)

	# TODO: find a way to eliminate the for-loops!
	partitionedMolecules = np.zeros((nObjects, nProcesses), np.bool)
	
	for moleculeIndex in np.arange(uniqueEntriesStructured.size):
		indexes = np.where(moleculeIndex == mapping)[0]
		randStream.numpyShuffle(indexes)

		for processIndex in np.arange(nProcesses):
			selectedIndexes = indexes[indexingRanges[moleculeIndex, processIndex]:indexingRanges[moleculeIndex, processIndex+1]]

			partitionedMolecules[
				selectedIndexes,
				processIndex] = True
	
	return partitionedMolecules

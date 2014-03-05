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

import numpy as np
import tables

import wholecell.utils.linear_programming as lp

class UniqueObjectsContainer(object):
	# TODO: move to separate file
	# TODO: generalize names (molecule -> object)
	'''
	UniqueObjectsContainer

	Essentially a dict of structured arrays, where the fields are attributes
	of the unique molecules.  Used for the unique molecules state and
	partitions.
	'''

	_defaultContainerAttributes = {
		'_isActive':'bool', # whether the row is an active entry
		'_wasDeleted':'bool', # whether the row was deleted in the last step
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

		self._savedAttributes = {} # objectName:list of attributes
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

			# Record which attributes are saved
			self._savedAttributes[objectName] = attributes.keys()
			self._savedAttributes[objectName].remove('_isActive') # only active molecules are saved, so this field is not needed

			# Give the tables accessible names
			self._tableNames[objectName] = objectName.replace(' ', '_')

		# TODO: alternate constructor for copying to partitions


	def moleculesNew(self, objectName, nMolecules, **attributes):
		arrayIndex = self._nameToArrayIndex[objectName]
		objectIndexes = self._getFreeIndexes(arrayIndex, nMolecules)

		array = self._arrays[arrayIndex]

		array['_isActive'][objectIndexes] = True

		for attrName, attrValue in attributes.viewitems():
			# NOTE: there is probably a non-loop solution to this, but the 'obvious' solution creates a copy instead of a view
			array[attrName][objectIndexes] = attrValue

		globalIndexes = self._getFreeIndexes(self._globalRefIndex, nMolecules)
		globalArray = self._arrays[self._globalRefIndex]
		globalArray['_isActive'][globalIndexes] = True
		globalArray['_arrayIndex'][globalIndexes] = arrayIndex
		globalArray['_objectIndex'][globalIndexes] = objectIndexes

		array['_globalIndex'][objectIndexes] = globalIndexes

		return self._molecules(arrayIndex, objectIndexes)


	def moleculeNew(self, objectName, **attributes):
		(molecule,) = self.moleculesNew(objectName, 1, **attributes) # NOTE: tuple unpacking

		return molecule


	def _getFreeIndexes(self, arrayIndex, nMolecules):
		freeIndexes = np.where(
			~self._arrays[arrayIndex]['_isActive']
			& ~self._arrays[arrayIndex]['_wasDeleted']
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


	def moleculesDel(self, molecules):
		for molecule in molecules:
			self.moleculeDel(molecule)


	def moleculeDel(self, molecule):
		self._clearEntry(molecule._arrayIndex, molecule._objectIndex)


	def _clearEntry(self, arrayIndex, objectIndex):
		array = self._arrays[arrayIndex]

		globalArray = self._arrays[self._globalRefIndex]
		globalArray[array[objectIndex]['_globalIndex']] = np.zeros(1, dtype = globalArray.dtype)
		globalArray[array[objectIndex]['_globalIndex']]['_wasDeleted'] = True

		array[objectIndex] = np.zeros(
			1,
			dtype = array.dtype
			)

		array[objectIndex]['_wasDeleted'] = True


	def _clearAll(self, arrayIndex):
		# NOTE: this a dangerous method, meant only to be called on load!
		self._arrays[arrayIndex] = np.zeros(0,
			dtype = self._arrays[arrayIndex].dtype)


	def evaluateQuery(self, objectName, **operations): # TODO: allow for queries over all or a subset of molecules
		arrayIndex = self._nameToArrayIndex[objectName]

		return self._molecules(
			arrayIndex,
			np.where(self._queryMolecules(arrayIndex, **operations))[0]
			)
	

	def _queryMolecules(self, arrayIndex, **operations):
		operations['_isActive'] = ('==', True)
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
				self._queryMolecules(query._arrayIndex, **query._operations)
				)[0]


	def queryNew(self, objectName, **operations):
		arrayIndex = self._nameToArrayIndex[objectName]

		query = _Query(self, arrayIndex, **operations)

		self._queries.append(query)

		# TODO: return an old query object if the query is the same

		return query


	def molecules(self, objectName):
		return self._molecules(self._nameToArrayIndex[objectName])


	def _molecules(self, arrayIndex, objectIndexes = None):
		return set(self._iterMolecules(arrayIndex, objectIndexes)) # TODO: return a set-like object that creates the _UniqueObject instances as needed


	def _moleculesByGlobalIndex(self, globalIndexes): # TODO: make global index ref the internal standard behavior?
		globalArray = self._arrays[self._globalRefIndex]

		return set(
			_UniqueObject(self, globalArray['_arrayIndex'][index], globalArray['_objectIndex'][index])
			for index in globalIndexes
			)


	def iterMolecules(self, objectName):
		return self._iterMolecules(self._nameToArrayIndex[objectName])


	def _iterMolecules(self, arrayIndex, objectIndexes = None):
		if objectIndexes is None:
			objectIndexes = np.where(self._arrays[arrayIndex]['_isActive'])[0]

		return (_UniqueObject(self, arrayIndex, objectIndex) for objectIndex in objectIndexes)



	def _timeIs(self, time):
		for array in self._arrays:
			array['_time'] = time


	def _flushDeleted(self):
		for array in self._arrays:
			array['_wasDeleted'] = False


	# TODO: fix saving/loading...
	# currently problematic because
	#	indexes won't line up on load which makes testing hard
	#	_wasDeleted property isnt saved
	#	global ref indexes point to the wrong spots

	# def pytablesCreate(self, h5file):
	# 	for arrayIndex, array in enumerate(self._arrays):
	# 		h5file.create_table(
	# 			h5file.root,
	# 			self._tableNames[arrayIndex],
	# 			array[self._savedAttributes[arrayIndex]].dtype,
	# 			title = self._objectNames[arrayIndex],
	# 			filters = tables.Filters(complevel = 9, complib = 'zlib')
	# 			)


	# def pytablesAppend(self, h5file, time):
	# 	for arrayIndex, array in enumerate(self._arrays):
	# 		table= h5file.get_node('/', self._tableNames[arrayIndex])

	# 		entries = array[self._queryMolecules(self._objectNames[array])][self._savedAttributes[arrayIndex]]

	# 		table.append(entries)

	# 		table.flush()


	# def pytablesLoad(self, h5file, timePoint):
	# 	for arrayIndex, tableName in enumerate(self._tableNames):
	# 		table = h5file.get_node('/', tableName)

	# 		entries = table[table[:]['_time'] == timePoint]

	# 		self._clearAll(arrayIndex)
	# 		indexes = self._getFreeIndexes(arrayIndex, entries.size)

	# 		for attrName in savedAttributes:
	# 			self._arrays[attrName][indexes] = entries[attrName]


	# TODO: compute mass
	# TODO: move to new file
	# TODO: implement molecule transfer between containers
	# TODO: __eq__ comparison for tests


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


	def molecules(self):
		return self._container._molecules(self._arrayIndex, self._objectIndexes)


	def iterMolecules(self):
		return self._container._iterMolecules(self._arrayIndex, self._objectIndexes)


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
		
		if not entry['_isActive']:
			raise Exception('Attempted to access an inactive molecule.')

		return entry[attribute]


	def attrIs(self, attribute, value):
		entry = self._container._arrays[self._arrayIndex][self._objectIndex]
		
		if not entry['_isActive']:
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

	nMolecules = objectRequestsArray.shape[0]
	nRequests = requestNumberVector.size
	nProcesses = requestProcessArray.shape[1]

	# Make into structured array to condense the problem into unique rows
	objectRequestsStructured = objectRequestsArray.view(
		dtype = [('', np.bool)] * nRequests)

	uniqueEntriesStructured, mapping = np.unique(objectRequestsStructured,
		return_inverse = True)

	uniqueEntries = uniqueEntriesStructured.view((np.bool, 
		len(uniqueEntriesStructured.dtype.names)))

	counts = np.bincount(mapping) # the number of each condensed molecule type

	nMoleculeTypes = counts.size

	# Some index mapping voodoo
	where0, where1 = np.where(uniqueEntries)

	nConnections = where0.size

	argsort = np.argsort(where1)

	moleculeToRequestConnections = np.zeros((nMoleculeTypes + nRequests,
		nConnections), np.int)

	upperIndices = (where0, np.arange(where1.size)[argsort])
	lowerIndices = (nMoleculeTypes + where1[argsort], np.arange(where1.size))
	# End voodoo

	moleculeToRequestConnections[upperIndices] = -1
	moleculeToRequestConnections[lowerIndices] = 1

	# Create the matrix and fill in the values
	matrix = np.zeros(
		(nMoleculeTypes + nRequests + nProcesses,
			nMoleculeTypes + nConnections + 2*nProcesses),
		np.int
		)

	# Molecule "boundary fluxes"
	matrix[:nMoleculeTypes, :nMoleculeTypes] = np.identity(nMoleculeTypes)

	# Flow from molecule type to request
	matrix[:nMoleculeTypes + nRequests,
		nMoleculeTypes:nMoleculeTypes+nConnections] = moleculeToRequestConnections

	# Flow from request to process
	matrix[nMoleculeTypes:nMoleculeTypes+nRequests,
		nMoleculeTypes+nConnections:nMoleculeTypes+nConnections+nProcesses][np.where(requestProcessArray)] = -requestNumberVector

	matrix[nMoleculeTypes + nRequests:,
		nMoleculeTypes+nConnections:nMoleculeTypes+nConnections+nProcesses] = np.identity(nProcesses)

	# Process "boundary fluxes"
	matrix[nMoleculeTypes + nRequests:,
		-nProcesses:] = -np.identity(nProcesses)

	# Create other linear programming parameters

	objective = np.zeros(matrix.shape[1], np.float)
	objective[-nProcesses:] = 1 # objective is to maximize process satisfaction
	# TODO: experiment with non-unity process weightings

	b = np.zeros(matrix.shape[0], np.float) # conservation law, i.e. b = 0 = Ax

	lowerBound = np.zeros(matrix.shape[1], np.float) # matrix is defined such that all values are >= 0

	upperBound = np.empty(matrix.shape[1], np.float)
	upperBound[:] = np.inf
	upperBound[:nMoleculeTypes] = counts # can use up to the total number of molecules
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

	unfixedCounts = -moleculeToRequestConnections[:nMoleculeTypes, :] * solution[nMoleculeTypes:nMoleculeTypes+nConnections]

	flooredCounts = np.floor(unfixedCounts) # Round down to prevent oversampling

	flooredProcessCounts = np.dot(
		np.dot(flooredCounts, moleculeToRequestConnections[nMoleculeTypes:, :].T),
		requestProcessArray
		)

	indexingRanges = np.c_[np.zeros(nMoleculeTypes), np.cumsum(flooredProcessCounts, 1)].astype(np.int)

	# TODO: find a way to eliminate the for-loops!
	partitionedMolecules = np.zeros((nMolecules, nProcesses), np.bool)
	
	for moleculeIndex in np.arange(uniqueEntriesStructured.size):
		indexes = np.where(moleculeIndex == mapping)[0]
		randStream.numpyShuffle(indexes)

		for processIndex in np.arange(nProcesses):
			selectedIndexes = indexes[indexingRanges[moleculeIndex, processIndex]:indexingRanges[moleculeIndex, processIndex+1]]

			partitionedMolecules[
				selectedIndexes,
				processIndex] = True
	
	return partitionedMolecules

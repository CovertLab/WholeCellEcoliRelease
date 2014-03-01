'''
unique_molecules.py

The UniqueMolecules State handles the identity and dynamic properties of unique
molecules in the simulation.  The attribute names and data types are imported
from the knowledge base.

The UniqueMolecules State instantiates a UniqueMoleculesContainer object,
which creates and manages the structured arrays in memory.

The UniqueMoleculesContainer uses _UniqueObject objects to present a clean 
interface to a specific molecule's attributes.

The UniqueMoleculesContainer also uses _Query objects to store and periodically
update queries which return sets of _UniqueObject that refer to molecules 
satisfying the query.
'''

import numpy as np
import tables

import wholecell.states.state

MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		'boundToChromosome':'bool',
		'chromosomeLocation':'uint32'
		}
	}

'''
TODO

Create an extra 'molecule' called _Global Reference, with attributes
_objectName (eventually, _arrayIndex)
_moleculeIndex

plus the default attributes.  Any add/remove action will update these entries.
Regular entries will get a back-reference index to their global entry.

Uses:
-unified access for queries over multiple molecules
-easier to wrap for other classes (sequence bound molecules)
-easy to track all unique molecules for loading, testing
-better support for indexed array references instead of dicts


'''

class UniqueMoleculesContainer(object):
	# TODO: move to separate file
	# TODO: generalize names (molecule -> object)
	'''
	UniqueMoleculesContainer

	Essentially a dict of structured arrays, where the fields are attributes
	of the unique molecules.  Used for the unique molecules state and
	partitions.
	'''

	_defaultContainerAttributes = {
		'_isActive':'bool', # whether the row is an active entry
		'_wasDeleted':'bool', # whether the row was deleted in the last step
		'_time':'uint64', # current time (important for saving)
		# '_globalIndex':'uint32'
		# '_massDifference':'float64' # dynamic mass difference
		}

	# _defaultMolecules = {
	# 	'_globalReference':{


	# 	}
	# }

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

		self._objectNames = sorted(self._objectAttributes.keys())

		for objectName, attributes in self._objectAttributes.viewitems():
			# Add the attributes used internally
			attributes.update(self._defaultContainerAttributes)

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

		indexes = self._getFreeIndexes(arrayIndex, nMolecules)

		array = self._arrays[arrayIndex]

		array['_isActive'][indexes] = True

		for attrName, attrValue in attributes.viewitems():
			# NOTE: there is probably a non-loop solution to this, but the 'obvious' solution creates a copy instead of a view
			array[attrName][indexes] = attrValue

		return self._molecules(arrayIndex, indexes)


	def moleculeNew(self, objectName, **attributes):
		(molecule,) = self.moleculesNew(objectName, 1, **attributes) # NOTE: tuple unpacking

		return molecule


	def moleculesDel(self, molecules):
		for molecule in molecules:
			self.moleculeDel(molecule)


	def moleculeDel(self, molecule):
		self._clearEntry(molecule._arrayIndex, molecule._index)


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


	def _clearEntry(self, arrayIndex, index):
		# this will probably be replaced
		array = self._arrays[arrayIndex]

		array[index] = np.zeros(
			1,
			dtype = array.dtype
			)

		array[index]['_wasDeleted'] = True


	def _clearAll(self, objectName):
		# NOTE: this a dangerous method, meant only to be called on load!
		self._arrays[self._nameToArrayIndex[objectName]] = np.zeros(0,
			dtype = self._arrays[self._nameToArrayIndex[objectName]].dtype)


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
			query._indexes = np.where(
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


	def _molecules(self, arrayIndex, indexes = None):
		return set(self._iterMolecules(arrayIndex, indexes)) # TODO: return a set-like object that creates the _UniqueObject instances as needed


	def iterMolecules(self, objectName):
		return self._iterMolecules(self._nameToArrayIndex[objectName])


	def _iterMolecules(self, arrayIndex, indexes = None):
		if indexes is None:
			indexes = np.where(self._arrays[arrayIndex]['_isActive'])[0]

		return (_UniqueObject(self, arrayIndex, index) for index in indexes)



	def _timeIs(self, time):
		for array in self._arrays:
			array['_time'] = time


	def _flushDeleted(self):
		for array in self._arrays:
			array['_wasDeleted'] = False


	def pytablesCreate(self, h5file):
		for objectName, savedAttributes in self._savedAttributes.viewitems():
			h5file.create_table(
				h5file.root,
				self._tableNames[objectName],
				self._nameToArray[objectName][savedAttributes].dtype,
				title = objectName,
				filters = tables.Filters(complevel = 9, complib = 'zlib')
				)


	def pytablesAppend(self, h5file, time):
		for objectName, savedAttributes in self._savedAttributes.viewitems():
			t = h5file.get_node('/', self._tableNames[objectName])

			entries = self._nameToArray[objectName][self._queryMolecules(objectName)][savedAttributes]

			t.append(entries)

			t.flush()


	def pytablesLoad(self, h5file, timePoint):
		for objectName, savedAttributes in self._savedAttributes.viewitems():
			t = h5file.get_node('/', self._tableNames[objectName])

			entries = t[t[:]['_time'] == timePoint]

			self._clearAll(objectName)
			indexes = self._getFreeIndexes(objectName, entries.size)

			for attrName in savedAttributes:
				self._nameToArray[objectName][attrName][indexes] = entries[attrName]


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

	__slots__ = ('_container', '_arrayIndex', '_operations', '_indexes')


	def __init__(self, container, arrayIndex, **operations):
		self._container = container
		self._arrayIndex = arrayIndex
		self._operations = operations
		self._indexes = None


	def molecules(self):
		return self._container._molecules(self._arrayIndex, self._indexes)


	def iterMolecules(self):
		return self._container._iterMolecules(self._arrayIndex, self._indexes)


	# TODO: sampling functions?  i.e. get N molecules
	# TODO: subqueries?


class _UniqueObject(object):
	'''
	_UniqueObject

	A wrapper around a row in a container, refering to a specific unique object.
	Primarily used as a way to manipulate individual molecules with a python
	object-like interface.
	'''
	
	__slots__ = ('_container', '_arrayIndex', '_index')


	def __init__(self, container, arrayIndex, index):
		self._container = container
		self._arrayIndex = arrayIndex
		self._index = index


	def attr(self, attribute):
		entry = self._container._arrays[self._arrayIndex][self._index]
		
		if not entry['_isActive']:
			raise Exception('Attempted to access an inactive molecule.')

		return entry[attribute]


	def attrIs(self, attribute, value):
		entry = self._container._arrays[self._arrayIndex][self._index]
		
		if not entry['_isActive']:
			raise Exception('Attempted to access an inactive molecule.')

		entry[attribute] = value


	def __hash__(self):
		return hash((self._container, self._arrayIndex, self._index))


	def __eq__(self, other):
		assert self._container == other._container, 'Molecule comparisons across UniqueMoleculesContainer objects not supported.'
		return self._arrayIndex == other._arrayIndex and self._index == other._index


	def __repr__(self):
		return '{}(..., {}, {})'.format(type(self).__name__, self._objectName, self._index)

	# TODO: method to get name


class UniqueMolecules(wholecell.states.state.State):
	'''
	UniqueMolecules

	State that tracks unique instances of molecules in the simulation, which 
	can have special dynamic attributes.
	'''


	def __init__(self, *args, **kwargs):
		self.meta = {
			'id':'UniqueMolecules',
			'name':'Unique Molecules',
			'dynamics':[],
			'units':{}
			}

		self.time = None

		self._container = None

		super(UniqueMolecules, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(UniqueMolecules, self).initialize(sim, kb)

		self.time = sim.states['Time']

		# TODO: use the updated KB object to get these properties

		self._container = UniqueMoleculesContainer(MOLECULE_ATTRIBUTES)

	
	def calcInitialConditions(self):
		# TODO: create a generalized calcInitialConditions routine as method of
		# the Simulation class, or as a separate function like fitSimulation

		# Add some molecules for testing save/load
		self._container.moleculesNew(
			'RNA polymerase',
			20,
			boundToChromosome = True, # just some example parameters
			chromosomeLocation = 50
			)


	def partition(self):
		# Set the correct time for saving purposes
		self._container._timeIs(self.time.value)

		# Clear out any deleted entries to make room for new molecules
		self._container._flushDeleted()

		# Update queries prior to gathering requests
		self._container.updateQueries()

		# TODO: actually partition


	def pytablesCreate(self, h5file, expectedRows):
		self._container.pytablesCreate(h5file)


	def pytablesAppend(self, h5file):
		self._container.pytablesAppend(h5file, self.time.value)


	def pytablesLoad(self, h5file, timePoint):
		self._container.pytablesLoad(h5file, timePoint)
	
	# TODO: partitioning

# TODO: partitions
# molecules created in a partition should be noted so they can be given a new, 
# permanent reference in the state

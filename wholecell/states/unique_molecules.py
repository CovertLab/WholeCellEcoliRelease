'''unique_molecules.py'''

import numpy as np
import tables

import wholecell.states.state

MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		'boundToChromosome':'bool',
		'chromosomeLocation':'uint32'
		}
	}

FRACTION_EXTEND_ENTRIES = 0.1 # fractional rate to increase number of entries in the structured array

QUERY_OPERATIONS = {
	'>':np.greater,
	'>=':np.greater_equal,
	'<':np.less,
	'<=':np.less_equal,
	'==':np.equal,
	'!=':np.not_equal
	}


class UniqueMoleculesContainer(object):
	'''
	UniqueMoleculesContainer

	Essentially a dict of structured arrays, where the fields are attributes
	of the unique molecules.  Used for the unique molecules state and
	partitions.
	'''

	defaultContainerAttributes = {
		'_isActive':'bool', # whether the row is an active entry
		'_wasDeleted':'bool', # whether the row was deleted in the last step
		'_time':'uint64', # current time (important for saving)
		# '_massDifference':'float64' # dynamic mass difference
		}


	def __init__(self, moleculeAttributes):
		self._moleculeAttributes = {}
		self._moleculeArrays = {}
		self._savedAttributes = {}
		self._queries = []
		self._tableNames = {}

		self._moleculeAttributes.update(moleculeAttributes)

		for moleculeName, attributes in self._moleculeAttributes.viewitems():
			# Add the attributes used internally
			attributes.update(self.defaultContainerAttributes)

			# Create the structured array
			self._moleculeArrays[moleculeName] = np.zeros(
				0, # start out empty
				dtype = [
					(attrName, attrType)
					for attrName, attrType in attributes.viewitems()
					]
				)

			# Record which attributes are saved
			self._savedAttributes[moleculeName] = attributes.keys()
			self._savedAttributes[moleculeName].remove('_isActive') # only active molecules are saved, so this field is not needed

			# Give the tables accessible names
			self._tableNames[moleculeName] = moleculeName.replace(' ', '_')

		# TODO: alternate constructor for copying to partitions


	def moleculesNew(self, moleculeName, nMolecules, **moleculeAttributes):
		indexes = self._getFreeIndexes(moleculeName, nMolecules)

		array = self._moleculeArrays[moleculeName]

		array['_isActive'][indexes] = True

		for attrName, attrValue in moleculeAttributes.viewitems():
			# NOTE: there is probably a non-loop solution to this, but the 'obvious' solution creates a copy instead of a view
			array[attrName][indexes] = attrValue

		return self.molecules(moleculeName, indexes)


	def moleculeNew(self, moleculeName, **moleculeAttributes):
		(molecule,) = self.moleculesNew(moleculeName, 1, **moleculeAttributes) # NOTE: tuple unpacking

		return molecule


	def moleculesDel(self, molecules):
		for molecule in molecules:
			self.moleculeDel(molecule)


	def moleculeDel(self, molecule):
		self._clearEntry(molecule._moleculeName, molecule._index)


	def _getFreeIndexes(self, moleculeName, nMolecules):
		freeIndexes = np.where(
			~self._moleculeArrays[moleculeName]['_isActive']
			& ~self._moleculeArrays[moleculeName]['_wasDeleted']
			)[0]

		if freeIndexes.size < nMolecules:
			oldEntries = self._moleculeArrays[moleculeName]
			oldSize = oldEntries.size

			newSize = oldSize + max(int(oldSize * FRACTION_EXTEND_ENTRIES), nMolecules)

			self._moleculeArrays[moleculeName] = np.zeros(
				newSize,
				dtype = oldEntries.dtype
				)
			
			self._moleculeArrays[moleculeName][:oldSize] = oldEntries

			freeIndexes = np.concatenate((freeIndexes, np.arange(oldSize, newSize)))

		return freeIndexes[:nMolecules]


	def _clearEntry(self, moleculeName, index):
		# this will probably be replaced
		array = self._moleculeArrays[moleculeName]

		array[index] = np.zeros(
			1,
			dtype = array.dtype
			)

		array[index]['_wasDeleted'] = True


	def _clearAll(self, moleculeName):
		# NOTE: this a dangerous method, meant only to be called on load!
		self._moleculeArrays[moleculeName] = np.zeros(0,
			dtype = self._moleculeArrays[moleculeName].dtype)


	def evaluateQuery(self, moleculeName, **operations): # TODO: allow for queries over all or a subset of molecules
		return self.molecules(moleculeName, np.where(self._queryMolecules(moleculeName, **operations))[0])
	

	def _queryMolecules(self, moleculeName, **operations):
		operations['_isActive'] = ('==', True)
		return reduce(
			np.logical_and,
			(
				QUERY_OPERATIONS[operator](
					self._moleculeArrays[moleculeName][attrName],
					queryValue
					)
				for attrName, (operator, queryValue) in operations.viewitems()
			)
		)


	def updateQueries(self):
		for query in self._queries:
			query._indexes = np.where(
				self._queryMolecules(query._moleculeName, **query._operations)
				)[0]


	def queryNew(self, moleculeName, **operations):
		query = _Query(self, moleculeName, **operations)

		self._queries.append(query)

		# TODO: return an old query object if the query is the same

		return query


	def molecules(self, moleculeName, _indexes = None):
		# NOTE: "_indexes" optional argument is for internal use; processes should use queries
		if _indexes is None:
			_indexes = np.where(self._moleculeArrays[moleculeName]['_isActive'])[0]

		return {_Molecule(self, moleculeName, index) for index in _indexes} # TODO: return a set-like object that create the _Molecule instances as needed


	def iterMolecules(self, moleculeName, _indexes = None):
		# NOTE: "_indexes" optional argument is for internal use; processes should use queries
		if _indexes is None:
			_indexes = np.where(self._moleculeArrays[moleculeName]['_isActive'])[0]

		return (_Molecule(self, moleculeName, index) for index in _indexes)


	def _timeIs(self, time):
		for array in self._moleculeArrays.viewvalues():
			array['_time'] = time


	def _flushDeleted(self):
		for array in self._moleculeArrays.viewvalues():
			array['_wasDeleted'] = False


	def pytablesCreate(self, h5file):
		for moleculeName, savedAttributes in self._savedAttributes.viewitems():
			h5file.create_table(
				h5file.root,
				self._tableNames[moleculeName],
				self._moleculeArrays[moleculeName][savedAttributes].dtype,
				title = moleculeName,
				filters = tables.Filters(complevel = 9, complib = 'zlib')
				)


	def pytablesAppend(self, h5file, time):
		for moleculeName, savedAttributes in self._savedAttributes.viewitems():
			t = h5file.get_node('/', self._tableNames[moleculeName])

			entries = self._moleculeArrays[moleculeName][self._queryMolecules(moleculeName)][savedAttributes]

			t.append(entries)

			t.flush()


	def pytablesLoad(self, h5file, timePoint):
		for moleculeName, savedAttributes in self._savedAttributes.viewitems():
			t = h5file.get_node('/', self._tableNames[moleculeName])

			entries = t[t[:]['_time'] == timePoint]

			self._clearAll(moleculeName)
			indexes = self._getFreeIndexes(moleculeName, entries.size)

			for attrName in savedAttributes:
				self._moleculeArrays[moleculeName][attrName][indexes] = entries[attrName]


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

	__slots__ = ('_container', '_moleculeName', '_operations', '_indexes')


	def __init__(self, container, moleculeName, **operations):
		self._container = container
		self._moleculeName = moleculeName
		self._operations = operations
		self._indexes = None


	def molecules(self):
		return self._container.molecules(self._moleculeName, self._indexes)


	def iterMolecules(self):
		return self._container.iterMolecules(self._moleculeName, self._indexes)


	# TODO: sampling functions?  i.e. get N molecules
	# TODO: subqueries?


class _Molecule(object):
	'''
	_Molecule

	A wrapper around a row in a container, refering to a specific molecule.
	'''
	
	__slots__ = ('_container', '_moleculeName', '_index')


	def __init__(self, container, moleculeName, index):
		self._container = container # TODO: cache reference to accessed array?
		self._moleculeName = moleculeName
		self._index = index


	def attr(self, attribute):
		entry = self._container._moleculeArrays[self._moleculeName][self._index]
		
		if not entry['_isActive']:
			raise Exception('Attempted to access an inactive molecule.')

		return entry[attribute]


	def attrIs(self, attribute, value):
		entry = self._container._moleculeArrays[self._moleculeName][self._index]
		
		if not entry['_isActive']:
			raise Exception('Attempted to access an inactive molecule.')

		entry[attribute] = value


	def __hash__(self):
		return hash((self._container, self._moleculeName, self._index))


	def __eq__(self, other):
		assert self._container == other._container, 'Molecule comparisons across UniqueMoleculesContainer objects not supported.'
		return self._moleculeName == other._moleculeName and self._index == other._index


	def __repr__(self):
		return '{}(..., {}, {})'.format(type(self).__name__, self._moleculeName, self._index)


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

'''unique_molecules.py'''

import numpy as np
import tables

import wholecell.states.state as wcState

MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		'boundToChromosome':'bool',
		'chromosomeLocation':'uint32'
		}
	}

FRACTION_EXTEND_ENTRIES = 0.1 # fractional rate to increase number of entries in the structured array

DEFAULT_STATE_ATTRIBUTES = { # attributes for local use by the state
	'_partitionedByOtherState':'bool', # whether the molecule is partitioned by a different state
	'_time':'uint64', # current time (important for saving)
	# '_massDifference':'float64' # dynamic mass difference
	}

DEFAULT_CONTAINER_ATTRIBUTES = { # attributes for local use by the container
	'_isActive':'bool', # whether the row is an active entry
	}

QUERY_OPERATIONS = {
	'>':np.greater,
	'>=':np.greater_equal,
	'<':np.less,
	'<=':np.less_equal,
	'==':np.equal,
	'!=':np.not_equal
	}

PYTABLES_TYPES = {
	'bool':tables.BoolCol,
	'uint32':tables.UInt32Col,
	'uint64':tables.UInt64Col
	}

class UniqueMoleculesContainer(object):
	'''
	UniqueMoleculesContainer

	Essentially a wrapper around a structured array, where the fields are 
	attributes of the unique molecules.  Used for the unique molecules state
	and partitions.
	'''

	def __init__(self, moleculeName, attributes):
		attributes.update(DEFAULT_CONTAINER_ATTRIBUTES)

		self._moleculeName = moleculeName
		self._attributes = attributes
		self._queries = []

		# Create the structured array needed to store the entries
		self._molecules = np.zeros(
			0, # start out empty
			dtype = [
				(attrName, attrType)
				for attrName, attrType in attributes.viewitems()
				]
			)

		# Record which attributes are saved
		self._savedAttributes = attributes.keys()
		self._savedAttributes.remove('_isActive') # only active molecules are saved
		self._tableName = self._moleculeName.replace(' ', '_')

		# TODO: alternate constructor for copying to partitions


	def moleculesNew(self, nMolecules, **moleculeAttributes):
		indexes = self._getFreeIndexes(nMolecules)

		self._molecules['_isActive'][indexes] = True

		for attrName, attrValue in moleculeAttributes.viewitems():
			# NOTE: there is probably a non-loop solution to this, but the 'obvious' solution creates a copy instead of a view
			self._molecules[attrName][indexes] = attrValue

		return self.molecules(indexes)


	def moleculeNew(self, **moleculeAttributes):
		(molecule,) = self.moleculesNew(1, **moleculeAttributes) # NOTE: tuple unpacking

		return molecule


	def moleculesDel(self, molecules):
		self._clearEntries([molecule._index for molecule in molecules])


	def moleculeDel(self, molecule):
		self._clearEntries(molecule._index)


	def _getFreeIndexes(self, nMolecules):
		freeIndexes = np.where(~self._molecules['_isActive'])[0]

		if freeIndexes.size < nMolecules:
			oldEntries = self._molecules
			oldSize = oldEntries.size

			newSize = oldSize + max(int(oldSize * FRACTION_EXTEND_ENTRIES), nMolecules)

			self._molecules = np.zeros(
				newSize,
				dtype = oldEntries.dtype
				)
			
			self._molecules[:oldSize] = oldEntries

			freeIndexes = np.concatenate((freeIndexes, np.arange(oldSize, newSize)))

		return freeIndexes[:nMolecules]


	def _clearEntries(self, indexes):
		# this will probably be replaced

		self._molecules[indexes] = np.zeros(
			1,
			dtype = self._molecules.dtype
			)


	def evaluateQuery(self, **operations):
		return self.molecules(np.where(self._queryMolecules(**operations))[0])
	

	def _queryMolecules(self, **operations):
		operations['_isActive'] = ('==', True)
		return reduce(
			np.logical_and,
			(
				QUERY_OPERATIONS[operator](
					self._molecules[attrName],
					queryValue
					)
				for attrName, (operator, queryValue) in operations.viewitems()
			)
		)


	def updateQueries(self):
		for query in self._queries:
			query._indexes = np.where(self._queryMolecules(**query._operations))[0]


	def queryNew(self, **operations):
		query = _Query(self, **operations)

		self._queries.append(query)

		return query


	def molecules(self, _indexes = None):
		# NOTE: "_indexes" optional argument is largely for internal use; processes should use queries
		if _indexes is None:
			_indexes = np.where(self._molecules['_isActive'])[0]

		return {_Molecule(self, index) for index in _indexes} # TODO: return a set-like object that create the _Molecule instances as needed


	def iterMolecules(self, _indexes = None):
		# NOTE: "_indexes" optional argument is largely for internal use; processes should use queries
		if _indexes is None:
			_indexes = np.where(self._molecules['_isActive'])[0]

		return (_Molecule(self, index) for index in _indexes) # TODO: return a set-like object that create the _Molecule instances as needed


	def pytablesCreate(self, h5file):
		# Create table
		h5file.create_table(
			h5file.root,
			self._tableName,
			self._molecules[self._savedAttributes].dtype,
			title = self._tableName,
			filters = tables.Filters(complevel = 9, complib = 'zlib')
			)


	def pytablesAppend(self, h5file, time):
		# Unfortunately this field *only* gets updated when appending, but 
		# it's not actually all that useful elsewhere
		self._molecules['_time'] = time

		t = h5file.get_node('/', self._tableName)

		entries = self._molecules[self._queryMolecules()][self._savedAttributes]

		t.append(entries)

		t.flush()


	def pytablesLoad(self, h5file, timePoint):
		raise NotImplementedError() # TODO


class _Query(object):
	'''
	_Query

	A reference to a query that can be updated by the parent container and 
	inspected by a process.
	'''

	__slots__ = ('_container', '_operations', '_indexes')

	def __init__(self, container, **operations):
		self._container = container
		self._operations = operations
		self._indexes = None


	def molecules(self):
		return self._container.molecules(self._indexes)

	def iterMolecules(self):
		return self._container.iterMolecules(self._indexes)

	# TODO: sampling functions?  i.e. get N molecules
	# TODO: subqueries?


class _Molecule(object):
	'''
	_Molecule

	A wrapper around a row in a container, refering to a specific molecule.
	'''
	
	__slots__ = ('_container', '_index')

	def __init__(self, container, index):
		self._container = container
		self._index = index


	def attr(self, attribute):
		return self._container._molecules[self._index][attribute]


	def attrIs(self, attribute, value):
		self._container._molecules[self._index][attribute] = value

	def __hash__(self):
		return hash((self._container, self._index))

	def __eq__(self, other):
		assert self._container == other._container, 'Molecule comparisons across UniqueMoleculesContainer objects not supported.'
		return self._index == other._index

	def __repr__(self):
		return '{}(..., {})'.format(type(self).__name__, self._index)


class UniqueMolecules(wcState.State):
	'''
	UniqueMolecules

	State that tracks unique instances of molecules in the simulation, which 
	can have special dynamic attributes.
	'''

	def __init__(self, *args, **kwargs):
		self.meta = {
			'id':'UniqueMolecules',
			'name':'Unique Molecules',
			'dyanamics':[],
			'units':{}
			}

		self.time = None

		self._containers = None

		super(UniqueMolecules, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(UniqueMolecules, self).initialize(sim, kb)

		self.time = sim.states['Time']

		# TODO: use the updated KB object to get these properties

		# NOTE: this is just to prevent modifying the underlying dictionaries with the default attributes
		moleculeAttributes = {
			moleculeName:attributes.copy()
			for moleculeName, attributes in MOLECULE_ATTRIBUTES.viewitems()
			}

		for attributes in moleculeAttributes.viewvalues():
			attributes.update(DEFAULT_STATE_ATTRIBUTES)

		self._containers = {
			moleculeName:UniqueMoleculesContainer(moleculeName, attributes)
			for moleculeName, attributes in moleculeAttributes.viewitems() # add default attrs here?
			}

	
	def calcInitialConditions(self):
		# TODO: create a generalized calcInitialConditions routine as method of
		# the Simulation class, or as a separate function like fitSimulation

		# Create some RNA polymerases with dummy properties
		self._containers['RNA polymerase'].moleculesNew(
			20,
			boundToChromosome = True, # just some example parameters
			chromosomeLocation = 50
			)

		# Run assertions
		rnaPolyContainer = self._containers['RNA polymerase']

		# Check the number of active entries
		activeEntries = rnaPolyContainer._molecules['_isActive']
		assert activeEntries.sum() == 20

		# Check that the active entries have the correct attribute value
		assert (rnaPolyContainer._molecules['boundToChromosome'][activeEntries] == True).all()
		assert (rnaPolyContainer._molecules['chromosomeLocation'][activeEntries] == 50).all()

		# Remove a few polymerases
		rnaPolyContainer._clearEntries(np.arange(5))

		# Check that the number of entries has decreased
		activeEntries = rnaPolyContainer._molecules['_isActive']
		assert activeEntries.sum() == 20 - 5

		# Raise a query
		active = rnaPolyContainer.evaluateQuery()

		# Check the query output
		assert len(active) == 15
		
		# Add a query
		boundToChromosome = rnaPolyContainer.queryNew(boundToChromosome = ('==', True))

		# Evaluate queries and test
		rnaPolyContainer.updateQueries()

		assert len(boundToChromosome.molecules()) == 15

		# Modify, update query, and assert
		rnaPolyContainer._clearEntries(np.arange(10))
		rnaPolyContainer.updateQueries()

		assert len(boundToChromosome.molecules()) == 10	

		# Check individual molecules
		for molecule in boundToChromosome.molecules():
			assert molecule.attr('boundToChromosome')

		# Set the attribute
		for molecule in boundToChromosome.molecules():
			molecule.attrIs('boundToChromosome', False)

		# Update the query and assert the changes
		rnaPolyContainer.updateQueries()

		assert not boundToChromosome.molecules() # should be empty

		# Clear out the molecules
		rnaPolyContainer.moleculesDel(rnaPolyContainer.molecules())

		assert not rnaPolyContainer.molecules()

		print 'All assertions passed.'

		# Add more molecules, for saving
		rnaPolyContainer.moleculesNew(
			20,
			boundToChromosome = True, # just some example parameters
			chromosomeLocation = 50
			)


	def pytablesCreate(self, h5file, expectedRows):
		for container in self._containers.viewvalues():
			container.pytablesCreate(h5file)


	def pytablesAppend(self, h5file):
		for container in self._containers.viewvalues():
			container.pytablesAppend(h5file, self.time.value)


	def pytablesLoad(self, h5file, timePoint):
		for container in self._containers.viewvalues():
			container.pytablesLoad(h5file, timePoint)
	
	# TODO: partitioning

# TODO: partitions
# molecules created in a partition should be noted so they can be given a new, 
# permanent reference in the state

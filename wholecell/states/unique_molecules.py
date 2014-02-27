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

DEFAULT_ATTRIBUTES = { # attributes for local use by the state
	'_isActive':'bool', # whether the row is an active entry
	'_partitionedByOtherState':'bool', # whether the molecule is partitioned by a different state
	# '_massDifference':'float64' # dynamic mass difference
	}

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

	Essentially a wrapper around a structured array, where the fields are 
	attributes of the unique molecules.  Used for the unique molecules state
	and partitions.
	'''

	def __init__(self, attributes):
		self._attributes = attributes
		self._queries = []

		# Create the structured array needed to store the entries
		self._molecules = np.zeros(
			0, # start out empty
			dtype = [
				(attributeName, attributeType)
				for attributeName, attributeType in attributes.viewitems()
				]
			)

		# TODO: alternate constructor for copying to partitions


	def moleculesNew(self, nMolecules, **moleculeAttributes):
		indexes = self._getFreeIndexes(nMolecules)

		self._molecules['_isActive'][indexes] = True

		for attribute, attrValue in moleculeAttributes.viewitems():
			# NOTE: there is probably a non-loop solution to this, but the 'obvious' solution creates a copy instead of a view
			self._molecules[attribute][indexes] = attrValue


	def moleculeNew(self, **moleculeAttributes):
		self.moleculesNew(1, **moleculeAttributes)


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
		operations['_isActive'] = ('==', True)
		return reduce(
			np.logical_and,
			(
				QUERY_OPERATIONS[operator](
					self._molecules[attribute],
					queryValue
					)
				for attribute, (operator, queryValue) in operations.viewitems()
			)
		)


	def updateQueries(self):
		for query in self._queries:
			query._indexes = np.where(self.evaluateQuery(**query._operations))[0]


	def queryNew(self, **operations):
		query = _Query(self, **operations)

		self._queries.append(query)

		return query


	def molecules(self, indexes = None):
		if indexes is None:
			indexes = np.where(self._molecules['_isActive'])[0]

		return {_Molecule(self, index) for index in indexes}


	# TODO: pytable create/save/load


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
			attributes.update(DEFAULT_ATTRIBUTES)

		self._containers = {
			moleculeName:UniqueMoleculesContainer(attributes)
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

		# Check the number of active entries
		activeEntries = self._containers['RNA polymerase']._molecules['_isActive']
		assert activeEntries.sum() == 20

		# Check that the active entries have the correct attribute value
		assert (self._containers['RNA polymerase']._molecules['boundToChromosome'][activeEntries] == True).all()
		assert (self._containers['RNA polymerase']._molecules['chromosomeLocation'][activeEntries] == 50).all()

		# Remove a few polymerases
		self._containers['RNA polymerase']._clearEntries(np.arange(5))

		# Check that the number of entries has decreased
		activeEntries = self._containers['RNA polymerase']._molecules['_isActive']
		assert activeEntries.sum() == 20 - 5

		# Raise a query
		active = self._containers['RNA polymerase'].evaluateQuery()
		boundToChromosome = self._containers['RNA polymerase'].evaluateQuery(boundToChromosome = ('==', True))
		multipleConditions = self._containers['RNA polymerase'].evaluateQuery(
			boundToChromosome = ('==', True),
			chromosomeLocation = ('>', 0)
			)
		notTrue = self._containers['RNA polymerase'].evaluateQuery(
			boundToChromosome = ('==', False),
			chromosomeLocation = ('>', 0)
			)

		# Check the query output
		assert active.sum() == 15
		assert (active == boundToChromosome).all()
		assert (active == multipleConditions).all()
		assert notTrue.sum() == 0

		# Add a query
		boundToChromosome = self._containers['RNA polymerase'].queryNew(boundToChromosome = ('==', True))

		# Evaluate queries and test
		self._containers['RNA polymerase'].updateQueries()

		assert len(boundToChromosome.molecules()) == 15

		# Modify, update query, and assert
		self._containers['RNA polymerase']._clearEntries(np.arange(10))
		self._containers['RNA polymerase'].updateQueries()

		assert len(boundToChromosome.molecules()) == 10	

		# Check individual molecules
		for molecule in boundToChromosome.molecules():
			assert molecule.attr('boundToChromosome')

		# Set the attribute
		for molecule in boundToChromosome.molecules():
			molecule.attrIs('boundToChromosome', False)

		# Update the query and assert the changes
		self._containers['RNA polymerase'].updateQueries()

		assert not boundToChromosome.molecules() # should be empty

		print 'All assertions passed.'
	
	# TODO: partitioning

# TODO: partitions
# molecules created in a partition should be noted so they can be given a new, 
# permanent reference in the state

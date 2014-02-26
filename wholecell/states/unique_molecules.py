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

N_ENTRIES = 10 # default number of entries in each structured array
FRACTION_EXTEND_ENTRIES = 0.1 # fractional rate to increase number of entries in the structured array

DEFAULT_ATTRIBUTES = [ # attributes for local use
	('_isActive', 'bool'), # whether the row is an active entry
	('_partitionedByOtherState', 'bool'), # whether the molecule is partitioned by a different state
	# ('_massDifference', 'float64') # dynamic mass difference
	]

QUERY_OPERATIONS = {
	'>':np.greater,
	'>=':np.greater_equal,
	'<':np.less,
	'<=':np.less_equal,
	'==':np.equal,
	'!=':np.not_equal
	}

class UniqueMolecules(wcState.State):
	def __init__(self, *args, **kwargs):
		self.meta = {
			'id':'UniqueMolecules',
			'name':'Unique Molecules',
			'dyanamics':[],
			'units':{}
			}

		self.time = None

		self._moleculeAttributes = {}
		self._uniqueMolecules = {}
		self._queries = []

		super(UniqueMolecules, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(UniqueMolecules, self).initialize(sim, kb)

		self.time = sim.states['Time']

		# TODO: use the updated KB object to get these properties

		# Collect attribute information from KB
		self._moleculeAttributes.update(MOLECULE_ATTRIBUTES)

		# Create the structured arrays needed to store the entries
		for moleculeName, attributes in self._moleculeAttributes.viewitems():
			self._uniqueMolecules[moleculeName] = np.zeros(
				N_ENTRIES,
				dtype = DEFAULT_ATTRIBUTES + [
					(attributeName, attributeType)
					for attributeName, attributeType in attributes.viewitems()
					]
				)

	
	def calcInitialConditions(self):
		# TODO: create a generalized calcInitialConditions routine as method of
		# the Simulation class, or as a separate function like fitSimulation

		# Create some RNA polymerases with dummy properties
		self.moleculesNew(
			'RNA polymerase', 20,
			boundToChromosome = True, # just some example parameters
			chromosomeLocation = 50
			)

		# Check the number of active entries
		activeEntries = self._uniqueMolecules['RNA polymerase']['_isActive']
		assert activeEntries.sum() == 20

		# Check that the active entries have the correct attribute value
		assert (self._uniqueMolecules['RNA polymerase']['boundToChromosome'][activeEntries] == True).all()
		assert (self._uniqueMolecules['RNA polymerase']['chromosomeLocation'][activeEntries] == 50).all()

		# Remove a few polymerases
		self._clearEntries('RNA polymerase', np.arange(5))

		# Check that the number of entries has decreased
		activeEntries = self._uniqueMolecules['RNA polymerase']['_isActive']
		assert activeEntries.sum() == 20 - 5

		# Raise a query
		active = self.query('RNA polymerase')
		boundToChromosome = self.query('RNA polymerase', boundToChromosome = ('==', True))
		multipleConditions = self.query(
			'RNA polymerase',
			boundToChromosome = ('==', True),
			chromosomeLocation = ('>', 0)
			)
		notTrue = self.query(
			'RNA polymerase',
			boundToChromosome = ('==', False),
			chromosomeLocation = ('>', 0)
			)

		# Check the query output
		assert (active == activeEntries).all()
		assert (active == boundToChromosome).all()
		assert (active == multipleConditions).all()
		assert (active[notTrue] == False).all()

		print 'All assertions passed.'


	def moleculeNew(self, moleculeName, **moleculeAttributes):
		self.moleculesNew(moleculeName, 1, **moleculeAttributes)


	def moleculesNew(self, moleculeName, nMolecules, **moleculeAttributes):
		indexes = self._getFreeIndexes(moleculeName, nMolecules)

		self._uniqueMolecules[moleculeName]['_isActive'][indexes] = True

		for attribute, attrValue in moleculeAttributes.viewitems():
			# NOTE: there is probably a non-loop solution to this, but the 'obvious' solution creates a copy instead of a view
			self._uniqueMolecules[moleculeName][attribute][indexes] = attrValue


	def _getFreeIndexes(self, moleculeName, nMolecules):
		freeIndexes = np.where(~self._uniqueMolecules[moleculeName]['_isActive'])[0]

		if freeIndexes.size < nMolecules:
			oldEntries = self._uniqueMolecules[moleculeName]
			oldSize = oldEntries.size

			newSize = oldSize + max(int(oldSize * FRACTION_EXTEND_ENTRIES), nMolecules)

			self._uniqueMolecules[moleculeName] = np.zeros(
				newSize,
				dtype = oldEntries.dtype
				)
			
			self._uniqueMolecules[moleculeName][:oldSize] = oldEntries

			freeIndexes = np.concatenate((freeIndexes, np.arange(oldSize, newSize)))

		return freeIndexes[:nMolecules]


	def _clearEntries(self, moleculeName, indexes):
		# this will probably be replaced

		self._uniqueMolecules[moleculeName][indexes] = np.zeros(
			1,
			dtype = self._uniqueMolecules[moleculeName].dtype
			)


	def query(self, moleculeName, **operations):
		operations['_isActive'] = ('==', True)
		return reduce(
			np.logical_and,
			(
				QUERY_OPERATIONS[operator](
					self._uniqueMolecules[moleculeName][attribute],
					queryValue
					)
				for attribute, (operator, queryValue) in operations.viewitems()
			)
			)

		# TODO: queries as objects?
		# TODO: return something more useful than a bool matrix

	# TODO: querying
	# TODO: partitioning
	# TODO: pytable create/save/load
	# TODO: accessors

# TODO: partitions
# challenges for partitions:
# * accessors for individual molecules, groups of molecules
# * merging new molecules

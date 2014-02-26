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

		# TODO

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

		self.moleculesNew(
			'RNA polymerase', 20,
			boundToChromosome = True, # just some example parameters
			chromosomeLocation = 50
			)

		activeEntries = self._uniqueMolecules['RNA polymerase']['_isActive']

		# Check the number of active entries
		assert activeEntries.sum() == 20

		# Check that the active entries have the correct attribute value
		assert (self._uniqueMolecules['RNA polymerase']['boundToChromosome'][activeEntries] == True).all()
		assert (self._uniqueMolecules['RNA polymerase']['chromosomeLocation'][activeEntries] == 50).all()


	def moleculeNew(self, moleculeName, **moleculeAttributes):
		self.moleculesNew(moleculeName, 1, **moleculeAttributes)


	def moleculesNew(self, moleculeName, nMolecules, **moleculeAttributes):
		indexes = self._getFreeIndexes(moleculeName, nMolecules)

		self._uniqueMolecules[moleculeName]['_isActive'][indexes] = True

		for attribute, attrValue in moleculeAttributes.viewitems():
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

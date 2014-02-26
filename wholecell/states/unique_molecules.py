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
				dtype = [('_isActive', 'bool'), ('partitionedByOtherState', 'bool')] + [
					(attributeName, attributeType)
					for attributeName, attributeType in attributes.viewitems()
					]
				)

			# _isActive: used internally to track which rows are active molecules
			# partitionedByOtherState: used internally and externally to track which molecules are partitioned by a different state

	
	def calcInitialConditions(self):
		# TODO: create a generalized calcInitialConditions routine as method of
		# the Simulation class, or as a separate function like fitSimulation

		self.moleculesNew(
			'RNA polymerase', 20,
			boundToChromosome = True,
			chromosomeLocation = 50
			)


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

			newSize = oldSize + min(int(oldSize * FRACTION_EXTEND_ENTRIES), nMolecules)

			self._uniqueMolecules[moleculeName] = np.zeros(
				newSize,
				dtype = oldEntries.dtype
				)
			
			self._uniqueMolecules[moleculeName][:oldSize] = oldEntries

			freeIndexes = np.concatenate((freeIndexes, np.arange(oldSize, newSize)))

		return freeIndexes[:nMolecules]

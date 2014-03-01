
import numpy as np
import tables

import wholecell.states.state
import wholecell.utils.unique_objects_container

N_BASES = 5000000 # TODO: from kb
MOLECULE_WIDTH = 50 # TODO: from kb

MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		}
	}

class SequenceBoundMolecules(object):
	# Special values for the state of the chromosome
	_inactive = 0
	_empty = 1

	_specialValues = np.array([_inactive, _empty])
	_offset = _specialValues.size

	_defaultMoleculesContainerAttributes = {
		'_sequenceBoundLocation':'uint32', # location of the molecule on the chromosome
		}

	
	def __init__(self):
		self._length = N_BASES
		
		self._array = np.empty(N_BASES, dtype = np.int32) # TODO: choose best dtype based on array size
		self._array[:] = self._inactive

		moleculeAttributes = {}
		for moleculeName, attributes in MOLECULE_ATTRIBUTES.viewitems():
			moleculeAttributes[moleculeName] = attributes.copy()
			moleculeAttributes[moleculeName].update(self._defaultMoleculesContainerAttributes)

		self._moleculesContainer =  wholecell.utils.unique_objects_container.UniqueObjectsContainer(
			moleculeAttributes)


	def moleculeNew(self, objectName, location, **attributes):
		molecule = self._moleculesContainer.moleculeNew(
			objectName,
			_sequenceBoundLocation = location,
			**attributes
			)

		self._setMoleculeLocation(molecule, location)

		return molecule


	def moleculeMove(self, molecule, newLocation):
		oldLocation = molecule.attr('_sequenceBoundLocation')
		self._array[oldLocation:oldLocation+MOLECULE_WIDTH] = self._empty

		self._setMoleculeLocation(molecule, newLocation)

		return molecule


	def _setMoleculeLocation(self, molecule, location):
		self._array[location:location+MOLECULE_WIDTH] = molecule.attr('_globalIndex') + self._offset
		molecule.attrIs('_sequenceBoundLocation', location)


	def moleculeDel(self, molecule):
		location = molecule.attr('_sequenceBoundLocation')
		self._array[location:location+MOLECULE_WIDTH] = self._empty

		self._moleculesContainer.moleculeDel(molecule)


	def molecules(self, start, stop):
		indexes = np.setdiff1d(self._array[start:stop], self._specialValues) - self._offset

		return self._moleculesContainer._moleculesByGlobalIndex(indexes)


	# TODO: figure out how molecule width info is going to be handled
	# TODO: circularly-permuted indexing
	# TODO: saving
	# TODO: update container time, flush deleted molecules, update queries?
	# TODO: handle/pass sequence, multiplicity
	# TODO: restrict binding to empty regions
	# TODO: methods for activating/inactivating regions
	# TODO: write tests


class Chromosome(wholecell.states.state.State):

	def __init__(self, *args, **kwargs):
		self.meta = {
			'id':'Chromosome',
			'name':'Chromosome',
			'dynamics':[],
			'units':{}
			}

		# self.time = None

		self._container = None

		super(Chromosome, self).__init__(*args, **kwargs)

	
	def initialize(self, sim, kb):
		super(Chromosome, self).initialize(sim, kb)

		self._container = SequenceBoundMolecules()


	def calcInitialConditions(self):
		# Randomly bind a bunch of RNA polymerases
		nBound = 0

		while nBound < 1000:
			bindAt = self.randStream.randi(N_BASES - MOLECULE_WIDTH)

			if not self._container.molecules(bindAt, bindAt + MOLECULE_WIDTH):
				self._container.moleculeNew('RNA polymerase', bindAt)

				nBound += 1


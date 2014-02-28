
import numpy as np
import tables

import wholecell.states.state
import wholecell.states.unique_molecules

N_BASES = 5000000 # TODO: from kb
MOLECULE_WIDTH = 50 # TODO: from kb

MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		'boundToChromosome':'bool',
		'chromosomeLocation':'uint32'
		}
	}

DEFAULT_SBM_ATTRIBUTES = {
	'_Location':'uint32',
	'_'
}

class SequenceBoundMolecules(object):
	# Special values for the state of the chromosome
	_inactive = 0
	_empty = 1

	_specialValues = np.array([_inactive, _empty])
	_offset = _specialValues.size

	_defaultMoleculesContainerAttributes = {
		'_sequenceBoundLocation':'uint32',
		'_sequenceBoundIndex':'uint32'
		}

	def __init__(self):
		# TODO: handle/pass sequence, multiplicity
		self._sequence = np.zeros(N_BASES, dtype = np.int8) # this really just needs to be a 2-bit integer...
		self._length = self._sequence.shape[0]
		
		self._array = np.empty(N_BASES, dtype = np.int32) # TODO: choose best dtype based on array size
		self._array[:] = self._empty

		moleculeAttributes = {}
		for moleculeName, attributes in MOLECULE_ATTRIBUTES:
			moleculeAttributes[moleculeName] = attributes.copy()
			moleculeAttributes[moleculeName].update(DEFAULT_SBM_ATTRIBUTES)

		self._moleculesContainer = UniqueMoleculesContainer(moleculeAttributes)
		self._boundMolecules = [] # TODO:
		# make this into a special UniqueMoleculeContainer-like object for holding references?
		# add a special molecule to the container? (I'm thinking this)


	def boundMolecules(self, start, stop):
		indexes = np.setdiff1d(self._array[start:stop], self._specialValues) - self._offset

		return {self._boundMolecules[index] for index in indexes}


	def boundMoleculeIs(self, molecule, start):
		index = self._getFreeIndexes(1)

		self._array[start:start+MOLECULE_WIDTH] = index
		molecule.attrIs('_sequenceBoundLocation', start)

		self._boundMolecules[index] = (molecule._moleculeName, molecule._index)
		molecule.attrIs('_sequenceBoundIndex', index)


	def _getFreeIndexes(self, n):
		raise NotImplementedError()



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
		super(UniqueMolecules, self).initialize(sim, kb)

		self._container = NucleotideBoundMolecules()


	def calcInitialConditions(self):
		pass


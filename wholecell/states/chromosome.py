
from __future__ import division

import numpy as np
import tables

import wholecell.states.state
import wholecell.utils.unique_objects_container

N_BASES = 5000000 # TODO: from kb
STRAND_MULTIPLICITY = 3
MOLECULE_WIDTH = 50 # TODO: from kb

MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		},
	'DNA polymerase':{
		},
	}

class SequenceBoundMolecules(object):
	# Special values for the state of the chromosome
	_inactive = 0
	_empty = 1
	_fork = 2 # the (inactive) location of a fork
	_end = 3 # the (inactive) location of the end of a strand, most likely where the fork is

	_specialValues = np.array([_inactive, _empty, _fork, _end])
	_offset = _specialValues.size

	_defaultMoleculesContainerAttributes = {
		'_sequenceBoundLocation':'uint32', # location of the molecule on the chromosome
		}

	_rootName = 'R' # The name of the root strand
	
	def __init__(self):
		self._length = N_BASES

		self._strandMultiplicity = STRAND_MULTIPLICITY
		self._nStrands = 1 + 2*(self._strandMultiplicity - 1)**2
		
		self._array = np.zeros((N_BASES, self._nStrands ), dtype = np.int32) # TODO: choose best dtype based on array size

		self._array[:, 0] = self._empty # Primary strand is always active

		moleculeAttributes = {}
		for moleculeName, attributes in MOLECULE_ATTRIBUTES.viewitems():
			moleculeAttributes[moleculeName] = attributes.copy()
			moleculeAttributes[moleculeName].update(self._defaultMoleculesContainerAttributes)

		self._moleculesContainer =  wholecell.utils.unique_objects_container.UniqueObjectsContainer(
			moleculeAttributes)

		self._buildStrandConnectivity()


	def _buildStrandConnectivity(self):
		strandNames = []
		strandNames.append(self._rootName)

		rootNames = strandNames
		for i in xrange(self._strandMultiplicity - 1):
			newNames = []

			for rootName in rootNames:
				newNames.append(rootName + 'A')
				newNames.append(rootName + 'B')

			strandNames.extend(newNames)
			rootNames = newNames

		self._strandNameToIndex = {name:ind for ind, name in enumerate(strandNames)}
		self._strandChildrenIndexes = [
			(self._strandNameToIndex[name + 'A'], self._strandNameToIndex[name + 'B'])
			if self._strandNameToIndex.has_key(name + 'A') else None
			for name in strandNames
			]

		self._strandParentIndex = [None]*len(strandNames)

		for parentIndex, childrenIndexes in enumerate(self._strandChildrenIndexes):
			if childrenIndexes is not None:
				self._strandParentIndex[childrenIndexes[0]] = parentIndex
				self._strandParentIndex[childrenIndexes[1]] = parentIndex


	def moleculeNew(self, objectName, location, **attributes):
		molecule = self._moleculesContainer.objectNew(
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

		self._moleculesContainer.objectDel(molecule)


	def molecules(self, start, stop):
		indexes = np.setdiff1d(self._array[start:stop], self._specialValues) - self._offset

		return self._moleculesContainer._objectsByGlobalIndex(indexes)


	def forkNew(self, strandName, start, stop):
		# TODO: assert fork does not exist
		# TODO: move things about randomly over region split or assert that it is empty

		strandParent = self._strandNameToIndex[strandName]
		try:
			strandChildA, strandChildB = self._strandChildrenIndexes[strandParent]

		except TypeError:
			raise Exception('No space allocated for strand {} to fork into'.format(strandName))


		self._array[start:stop, strandParent] = self._inactive
		self._array[start+1, strandParent] = self._fork
		self._array[stop-1, strandParent] = self._fork

		self._array[start:stop, strandChildA] = self._empty
		self._array[start-1, strandChildA] = self._end
		self._array[stop+1, strandChildA] = self._end

		self._array[start:stop, strandChildB] = self._empty
		self._array[start-1, strandChildB] = self._end
		self._array[stop+1, strandChildB] = self._end


	

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
		# Add replication fork

		origin = N_BASES // 2

		forkStart = origin - N_BASES // 4
		forkStop = origin + N_BASES // 4

		self._container.forkNew('R', forkStart, forkStop)

		
		import ipdb; ipdb.set_trace()



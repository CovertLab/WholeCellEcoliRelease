'''
unique_molecules.py

The UniqueMolecules State handles the identity and dynamic properties of unique
molecules in the simulation.  The attribute names and data types are imported
from the knowledge base.

The UniqueMolecules State instantiates a UniqueObjectsContainer object, which 
creates and manages the structured arrays in memory.
'''

import numpy as np
import tables

import wholecell.states.state
import wholecell.states.partition
import wholecell.utils.unique_objects_container


MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		'boundToChromosome':'bool',
		'chromosomeLocation':'uint32',
		'_partitionedProcess':'uint32' # TODO: assign to every as default instead of in this fake KB
		}
	}


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
		self.partitionClass = UniqueMoleculesPartition

		self._container = None

		super(UniqueMolecules, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(UniqueMolecules, self).initialize(sim, kb)

		self.time = sim.states['Time']

		# TODO: use the updated KB object to get these properties

		self._container = wholecell.utils.unique_objects_container.UniqueObjectsContainer(
			MOLECULE_ATTRIBUTES)

		for iProcess, partition in enumerate(self.partitions.viewvalues()):
			partition._processIndexIs(iProcess + 1)

	
	def calcInitialConditions(self):
		# TODO: create a generalized calcInitialConditions routine as method of
		# the Simulation class, or as a separate function like fitSimulation

		# # Add some molecules for testing save/load
		# self._container.moleculesNew(
		# 	'RNA polymerase',
		# 	20,
		# 	boundToChromosome = True, # just some example parameters
		# 	chromosomeLocation = 50
		# 	)

		pass


	def partition(self):
		# Set the correct time for saving purposes
		self._container._timeIs(self.time.value)

		# Clear out any deleted entries to make room for new molecules
		self._container._flushDeleted()

		# Update queries prior to gathering requests
		self._container.updateQueries()

		# Gather requests
		arrays = []
		counts = []
		processReference = []

		for iProcess, partition in enumerate(self.partitions.viewvalues()):
			partitionArrays, partitionCounts = partition.request()

			nRequests = len(partitionArrays)

			arrays.extend(partitionArrays)
			counts.extend(partitionCounts)
			processReference.extend([iProcess]*nRequests)

		nTotalRequests = len(processReference)
		nProcesses = len(self.partitions)

		# Format requests into appropriate parameters
		objectRequestsArray = np.vstack(arrays).transpose().copy() # must copy to fix memory order
		requestNumberVector = np.array(counts)
		requestProcessArray = (np.tile(np.arange(nProcesses), (nTotalRequests, 1)).T
			== np.array(processReference)).transpose()

		# import ipdb; ipdb.set_trace()

		partitionedMolecules = wholecell.utils.unique_objects_container._partition(
			objectRequestsArray, requestNumberVector, requestProcessArray, self.randStream)

		for iProcess, partition in enumerate(self.partitions.viewvalues()):
			molecules = self._container._objectsByGlobalIndex(
				np.where(partitionedMolecules[:, iProcess])[0]
				)

			for molecule in molecules:
				molecule.attrIs('_partitionedProcess', iProcess+1)


	def queryNew(self, moleculeName, **operations):
		return self._container.queryNew(moleculeName, **operations)


	def pytablesCreate(self, h5file, expectedRows):
		# self._container.pytablesCreate(h5file)
		pass


	def pytablesAppend(self, h5file):
		# self._container.pytablesAppend(h5file, self.time.value)
		pass


	def pytablesLoad(self, h5file, timePoint):
		# self._container.pytablesLoad(h5file, timePoint)
		pass
	


class UniqueMoleculesPartition(wholecell.states.partition.Partition):
	def __init__(self, *args, **kwargs):
		self._requestsEmpty()

		self._processIndex = None

		super(UniqueMoleculesPartition, self).__init__(*args, **kwargs)


	def _processIndexIs(self, value):
		self._processIndex = value


	def _requestsEmpty(self):
		self.requestArrays = []
		self.requestCounts = []


	def request(self):
		self._requestsEmpty()
		self._process.requestUniqueMolecules()

		return self.requestArrays, self.requestCounts


	def requestByMolecules(self, nMolecules, molecules):
		if nMolecules > 0:
			self.requestCounts.append(nMolecules)

			container = self._state._container

			size = container._arrays[container._globalRefIndex].size

			array = np.zeros(size, np.bool)

			indexes = np.array([molecule.attr('_globalIndex') for molecule in molecules])

			array[indexes] = True

			self.requestArrays.append(array)


	def evaluateQuery(self, moleculeName, **operations):
		operations['_partitionedProcess'] = ('==', self._processIndex)

		return self._state._container.evaluateQuery(moleculeName, **operations)


	def moleculesDel(self, molecules):
		self._state._container.objectsDel(molecules)


	def moleculeDel(self, molecule):
		self._state._container.objectDel(molecule)


	def moleculesNew(self, moleculeName, nMolecules, **attributes):
		attributes['_partitionedProcess'] = self._processIndex

		self._state._container.objectsNew(moleculeName, nMolecules, **attributes)


	def moleculeNew(self, moleculeName, **attributes):
		attributes['_partitionedProcess'] = self._processIndex

		self._state._container.objectNew(moleculeName, **attributes)




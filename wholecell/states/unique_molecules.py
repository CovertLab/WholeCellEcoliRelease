'''
unique_molecules.py

The UniqueMolecules State handles the identity and dynamic properties of unique
molecules in the simulation.  The attribute names and data types are imported
from the knowledge base.

The UniqueMolecules State instantiates a UniqueObjectsContainer object, which 
creates and manages the structured arrays in memory.
'''

from __future__ import division

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

		self._container = None

		super(UniqueMolecules, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(UniqueMolecules, self).initialize(sim, kb)

		self.time = sim.states['Time']

		# TODO: use the updated KB object to get these properties

		self._container = wholecell.utils.unique_objects_container.UniqueObjectsContainer(
			MOLECULE_ATTRIBUTES)

	
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


	def updateQueries(self):
		self._container.updateQueries()

		super(UniqueMolecules, self).updateQueries()


	def partition(self):
		# Set the correct time for saving purposes
		self._container._timeIs(self.time.value)

		# Clear out any deleted entries to make room for new molecules
		self._container._flushDeleted()
		
		# Gather requests
		nMolecules = self._container._arrays[self._container._globalRefIndex].size
		nViews = len(self._views)

		objectRequestsArray = np.zeros((nMolecules, nViews), np.bool)
		requestNumberVector = np.zeros(nViews, np.int64)
		requestProcessArray = np.zeros((nViews, self._nProcesses), np.bool)

		for viewIndex, view in enumerate(self._views):
			objectRequestsArray[view._queryObject._globalIndexes(), viewIndex] = True

			requestNumberVector[viewIndex] = view._request()

			requestProcessArray[viewIndex, view._processIndex] = True

		partitionedMolecules = wholecell.utils.unique_objects_container._partition(
			objectRequestsArray, requestNumberVector, requestProcessArray, self.randStream)

		for view in self._views:
			molecules = self._container._objectsByGlobalIndex(
				np.where(partitionedMolecules[:, view._processIndex])[0]
				)

			for molecule in molecules:
				molecule.attrIs(
					'_partitionedProcess',
					view._processIndex + 1 # "0", being the default, is reserved for unpartitioned molecules
					)


	def pytablesCreate(self, h5file, expectedRows):
		# self._container.pytablesCreate(h5file)
		pass


	def pytablesAppend(self, h5file):
		# self._container.pytablesAppend(h5file, self.time.value)
		pass


	def pytablesLoad(self, h5file, timePoint):
		# self._container.pytablesLoad(h5file, timePoint)
		pass

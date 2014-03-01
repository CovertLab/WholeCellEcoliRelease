'''
unique_molecules.py

The UniqueMolecules State handles the identity and dynamic properties of unique
molecules in the simulation.  The attribute names and data types are imported
from the knowledge base.

The UniqueMolecules State instantiates a UniqueMoleculesContainer object,
which creates and manages the structured arrays in memory.

The UniqueMoleculesContainer uses _UniqueObject objects to present a clean 
interface to a specific molecule's attributes.

The UniqueMoleculesContainer also uses _Query objects to store and periodically
update queries which return sets of _UniqueObject that refer to molecules 
satisfying the query.
'''

import numpy as np
import tables

import wholecell.states.state
import wholecell.utils.unique_objects_container


MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		'boundToChromosome':'bool',
		'chromosomeLocation':'uint32'
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

		# Add some molecules for testing save/load
		self._container.moleculesNew(
			'RNA polymerase',
			20,
			boundToChromosome = True, # just some example parameters
			chromosomeLocation = 50
			)


	def partition(self):
		# Set the correct time for saving purposes
		self._container._timeIs(self.time.value)

		# Clear out any deleted entries to make room for new molecules
		self._container._flushDeleted()

		# Update queries prior to gathering requests
		self._container.updateQueries()

		# TODO: actually partition


	def pytablesCreate(self, h5file, expectedRows):
		# self._container.pytablesCreate(h5file)
		pass


	def pytablesAppend(self, h5file):
		# self._container.pytablesAppend(h5file, self.time.value)
		pass


	def pytablesLoad(self, h5file, timePoint):
		# self._container.pytablesLoad(h5file, timePoint)
		pass
	
	# TODO: partitioning

# TODO: partitions
# molecules created in a partition should be noted so they can be given a new, 
# permanent reference in the state

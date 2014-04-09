
from __future__ import division

import numpy as np
import tables

import wholecell.states.state
from wholecell.containers.chromosome_container import ChromosomeContainer, ChromosomeContainerException

N_BASES = 5000000 # TODO: from kb
STRAND_MULTIPLICITY = 3 # TODO: estimate somehow? from kb?

MOLECULE_ATTRIBUTES = {
	'RNA polymerase':{
		},
	'DNA polymerase':{
		},
	}


# TODO: design views, queries, requests, etc

class Chromosome(wholecell.states.state.State):

	def __init__(self, *args, **kwargs):
		self.meta = {
			'id':'Chromosome',
			'name':'Chromosome',
			'dynamics':[],
			'units':{}
			}

		self.container = None

		self._partitioningHierarchy = ['ToyReplication', 'ToyTranscription']
		self._processToHierarchy = {name:i
			for i, name in enumerate(self._partitioningHierarchy)}

		super(Chromosome, self).__init__(*args, **kwargs)

	
	def initialize(self, sim, kb, kb2):
		super(Chromosome, self).initialize(sim, kb, kb2)

		self.container = ChromosomeContainer(N_BASES, STRAND_MULTIPLICITY,
			MOLECULE_ATTRIBUTES)

		self.time = sim.states['Time']


	def calcInitialConditions(self):
		# Replicate some regions
		forksRoot = self.container.divideRegion(
			self.container.rootStrand(),
			N_BASES // 6,
			5 * N_BASES //6
			)

		forksA = self.container.divideRegion(
			self.container.rootStrand() + 'A',
			N_BASES // 3,
			2 * N_BASES //3
			)

		forksB = self.container.divideRegion(
			self.container.rootStrand() + 'B',
			N_BASES // 3,
			2 * N_BASES //3
			)

		# Bind DNA polymerase
		forks = forksRoot + forksA + forksB

		for fork in forks:
			dnaPoly = self.container.moleculeNew('DNA polymerase')

			self.container.moleculeLocationIsFork(dnaPoly, fork, 30, 20)

		# Bind RNA polymerases randomly to the root strand

		unboundMolecules = [self.container.moleculeNew('RNA polymerase')
			for i in xrange(100)]

		while unboundMolecules:
			position = self.randStream.randi(N_BASES)
			orientation = ['+', '-'][self.randStream.randi(1)]

			try:
				self.container.moleculeLocationIs(
					unboundMolecules[-1],
					self.container.rootStrand(),
					position,
					orientation,
					30,
					20
					)

			except ChromosomeContainerException:
				continue

			else:
				unboundMolecules.pop()


	def partition(self):
		# Set the correct time for saving purposes
		self.container.timeIs(self.time.value)

		# TODO: flush deleted

		# First implementation: hierarchical
		# - Processes' requests are allocated in a fixed order
		# - Regions overlapping with previously allocated nucleotides are invalid

		# Order views by process

		viewGroups = [[] for process in self._partitioningHierarchy]

		for view in self._views: # NOTE: this needn't be reestablished every time step
			process = view._processId
			processIndex = self._processToHierarchy[process]
			viewGroups[processIndex].append(view)

		allocatedRegions = []

		# this logic is stupid and i hate it
		for viewGroup in viewGroups:
			regions = []
			for view in viewGroup:
				for region in view.requestedRegions():
					indexes = region.indexes()
					for allocatedRegion in allocatedRegions:
						if np.lib.arraysetops.intersect1d(
								indexes,
								allocatedRegion.indexes()
								).any():

							break

					else:
						regions.append(region)

				view.allocateRegions(regions)

			allocatedRegions.extend(regions)


	def calculate(self):
		# print 'Chromosome multiplicity is ', (self.container._array != self.container._inactive).sum() / N_BASES
		pass


	def pytablesCreate(self, h5file, expectedRows):
		self.container.pytablesCreate(h5file)


	def pytablesAppend(self, h5file):
		self.container.pytablesAppend(h5file)


	def pytablesLoad(self, h5file, timePoint):
		self.container.pytablesLoad(h5file, timePoint)


"""
Complexation

Macromolecular complexation sub-model. Encodes molecular simulation of macromolecular complexation

TODO:
- allow for shuffling when appropriate (maybe in another process)
- handle protein complex dissociation

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/4/2013
"""

from __future__ import absolute_import, division, print_function

import numpy as np
from arrow import StochasticSystem

import wholecell.processes.process

# Maximum unsigned int value + 1 for randint() to seed srand from C stdlib
RAND_MAX = 2**31

class Complexation(wholecell.processes.process.Process):
	""" Complexation """

	_name = "Complexation"

	# Constructor
	def __init__(self):

		super(Complexation, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(Complexation, self).initialize(sim, sim_data)

		# Create matrices and vectors that describe reaction stoichiometries
		self.stoichMatrix = sim_data.process.complexation.stoich_matrix().astype(np.int64)

		# semi-quantitative rate constants
		self.rates = sim_data.process.complexation.rates

		# build stochastic system simulation
		seed = self.randomState.randint(RAND_MAX)
		self.system = StochasticSystem(self.stoichMatrix.T, random_seed=seed)

		# Build views
		moleculeNames = sim_data.process.complexation.molecule_names
		self.molecules = self.bulkMoleculesView(moleculeNames)


	def calculateRequest(self):
		moleculeCounts = self.molecules.total_counts()

		result = self.system.evolve(
			self._sim.timeStepSec(), moleculeCounts, self.rates)
		updatedMoleculeCounts = result['outcome']

		self.molecules.requestIs(np.fmax(moleculeCounts - updatedMoleculeCounts, 0))


	def evolveState(self):
		moleculeCounts = self.molecules.counts()

		result = self.system.evolve(
			self._sim.timeStepSec(), moleculeCounts, self.rates)
		updatedMoleculeCounts = result['outcome']
		events = result['occurrences']

		self.molecules.countsIs(updatedMoleculeCounts)

		# Write outputs to listeners
		self.writeToListener("ComplexationListener", "complexationEvents", events)

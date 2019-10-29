#!/usr/bin/env python

"""
CellDivision

Cell division listener. Checks for cell division criteria..

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/18/2016
"""

# TODO: generalize this logic for use with a generic simulation

from __future__ import division

import numpy as np

import wholecell.listeners.listener
from wholecell.utils import units

class CellDivision(wholecell.listeners.listener.Listener):
	""" CellDivision """

	_name = 'CellDivision'

	# Constructor
	def __init__(self, *args, **kwargs):
		# References to other states
		self.internal_states = None

		# NOTE: molecule weight is converted to femtograms/molecule from
		# grams/mol in BulkMolecules
		self.massUnits = 'fg'

		super(CellDivision, self).__init__(*args, **kwargs)

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(CellDivision, self).initialize(sim, sim_data)

		self.internal_states = sim.internal_states

		self.waterIndex = sim_data.submassNameToIndex["water"]

		# Set total mass that should be added to cell
		# This is an approximation for length
		self.expectedDryMassIncreaseDict = sim_data.expectedDryMassIncreaseDict

		# Set initial values

		self.setInitial = False
		self.dryMass = 0.0
		# TODO: set initial masses based on some calculations of the expected
		# mother cell (divided by two) in the last time step

		self.massCoeff = 1.0
		if sim._massDistribution:
			self.massCoeff = sim.randomState.normal(loc = 1.0, scale = 0.1)

		# View on full chromosomes
		self.fullChromosomeView = self.internal_states['BulkMolecules'].container.countView('CHROM_FULL[c]')
		self.partialChromosomeView = self.internal_states['BulkMolecules'].container.countsView(self.internal_states['BulkMolecules'].divisionIds['partialChromosome'])
		self.fullChromosomeView = self.internal_states['BulkMolecules'].container.countView(sim_data.moleculeIds.fullChromosome)
		self.uniqueMoleculeContainer = self.internal_states['UniqueMolecules'].container

		self.divisionMassMultiplier = 1.
		if sim._massDistribution:
			self.divisionMassMultiplier = sim.randomState.normal(loc = 1.0, scale = 0.1)

		self.d_period_division = False
		if sim._dPeriodDivision:
			self.d_period_division = True

	def update(self):
		masses = sum(state.mass() for state in self.internal_states.itervalues())

		postEvolveMasses = masses[1, ...]

		self.cellMass = postEvolveMasses.sum() # sum over all dimensions
		submasses = postEvolveMasses.sum(axis = 0) # sum over the processes

		self.waterMass = submasses[self.waterIndex]
		self.dryMass = self.cellMass - self.waterMass

		if not self.setInitial:
			self.setInitial = True
			self.dryMassInitial = self.dryMass

		partial_chromosome_counts = self.partialChromosomeView.counts()
		uneven_counts = partial_chromosome_counts - partial_chromosome_counts.min()

		# Ends simulation once D period has occured after chromosome termination

		if self.d_period_division:
			fullChrom = self.uniqueMoleculeContainer.objectsInCollection("fullChromosome")
			if len(fullChrom):
				division_times = fullChrom.attr("division_time")
				divide_at_time = division_times.min()

				if self.time() >= divide_at_time and not uneven_counts.any():
					fullChrom.delByIndexes(np.where(division_times == divide_at_time)[0])
					self._sim.cellCycleComplete()
		else:
			# End simulation once the mass of an average cell is
			# added to current cell.
			current_nutrients = self._external_states['Environment'].nutrients
			if self.dryMass - self.dryMassInitial >= self.expectedDryMassIncreaseDict[current_nutrients].asNumber(units.fg) * self.divisionMassMultiplier:
				if not uneven_counts.any():
					self._sim.cellCycleComplete()

#!/usr/bin/env python

"""
CellDivision

Cell division listener. Checks for cell division criteria.

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

		# Get container for unique molecules
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

		# Ends simulation once D period has occurred after chromosome termination
		if self.d_period_division:
			# Get all existing full chromosomes
			full_chromosomes = self.uniqueMoleculeContainer.objectsInCollection("fullChromosome")

			# If there are two or more full chromosomes,
			if len(full_chromosomes) >= 2:
				# Extract attributes from existing full chromosomes
				division_time, chromosomeIndex = full_chromosomes.attrs(
					"division_time", "chromosomeIndex"
					)

				# Set division time to be the time set by the chromosome with
				# index 1
				divide_at_time = division_time[chromosomeIndex == 1][0]

				if self.time() >= divide_at_time:
					self._sim.cellCycleComplete()
		else:
			# End simulation once the mass of an average cell is
			# added to current cell.
			current_nutrients = self._external_states['Environment'].nutrients
			if self.dryMass - self.dryMassInitial >= self.expectedDryMassIncreaseDict[current_nutrients].asNumber(units.fg) * self.divisionMassMultiplier:
				self._sim.cellCycleComplete()

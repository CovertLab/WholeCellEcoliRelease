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
from wholecell.containers.unique_objects_container import Access

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

		# If D_PERIOD_DIVISION is set to True, the cell should divide after a
		# fixed length of time (D period) has passed after the completion of
		# chromosome replication. At completion of replication, we calculate
		# this division time by adding the D period length to the current time,
		# and storing this value as an attribute of the full_chromosome
		# molecule that gets formed at termination of replication. It is
		# possible for one more round of replication to be completed before
		# the cell divides if the cell is growing fast - so we take the minimum
		# division_time value from all existing full chromosomes, and mark the
		# chromosome if the chromosome has already induced division to avoid
		# double counting.
		if self.d_period_division:
			# Get all existing full chromosomes
			full_chromosomes = self.uniqueMoleculeContainer.objectsInCollection(
				"fullChromosome", access=[Access.EDIT])

			# If there are two or more full chromosomes,
			if len(full_chromosomes) >= 2:
				# Extract attributes from existing full chromosomes
				division_time, has_induced_division = full_chromosomes.attrs(
					"division_time", "has_induced_division"
					)

				# Set division time to be the minimum division_time of the
				# chromosome that have not yet induced cell division
				divide_at_time = division_time[~has_induced_division].min()

				if self.time() >= divide_at_time:
					self._sim.cellCycleComplete()

					# Mark the chromosome that induced division
					divide_at_time_index = np.where(
						division_time == divide_at_time)[0][0]
					has_induced_division[divide_at_time_index] = True
					full_chromosomes.attrIs(
						has_induced_division=has_induced_division
						)

		else:
			# End simulation once the mass of an average cell is
			# added to current cell.
			current_media_id = self._external_states['Environment'].current_media_id
			if self.dryMass - self.dryMassInitial >= self.expectedDryMassIncreaseDict[current_media_id].asNumber(units.fg) * self.divisionMassMultiplier:
				self._sim.cellCycleComplete()

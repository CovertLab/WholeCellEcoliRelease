"""
CellDivision process

- Flags the cell for division when a preset division criterion has been met

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/6/2020
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.processes.process
from wholecell.utils import units

class CellDivision(wholecell.processes.process.Process):
	""" CellDivision """

	_name = "CellDivision"

	# Constructor
	def __init__(self):
		super(CellDivision, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(CellDivision, self).initialize(sim, sim_data)

		self.d_period_division = sim._dPeriodDivision

		if self.d_period_division:
			self.full_chromosomes = self.uniqueMoleculesView('full_chromosome')

		else:
			# Flag for initial mass
			self.initial_mass_set = False

			# Set total mass that should be added to cell for division (this is
			# an approximation for length)
			self.expectedDryMassIncreaseDict = sim_data.expectedDryMassIncreaseDict
			self.internal_states = sim.internal_states
			self.water_index = sim_data.submassNameToIndex["water"]

			if sim._massDistribution:
				self.division_mass_multiplier = sim.randomState.normal(loc=1.0, scale=0.1)
			else:
				self.division_mass_multiplier = 1.


	def calculateRequest(self):
		if self.d_period_division:
			self.full_chromosomes.request_access(self.EDIT_ACCESS)


	def evolveState(self):
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
			if self.full_chromosomes.total_counts() >= 2:
				# Extract attributes from existing full chromosomes
				division_time, has_triggered_division = self.full_chromosomes.attrs(
					"division_time", "has_triggered_division"
					)

				# Set division time to be the minimum division_time of the
				# chromosome that has not yet triggered cell division
				divide_at_time = division_time[~has_triggered_division].min()

				if self.time() >= divide_at_time:
					self._sim.cellCycleComplete()

					# Mark the chromosome that induced division
					divide_at_time_index = np.where(
						division_time == divide_at_time)[0][0]
					has_triggered_division[divide_at_time_index] = True
					self.full_chromosomes.attrIs(
						has_triggered_division=has_triggered_division
						)

		else:
			# Calculate dry mass before timestep
			all_submasses_pre_timestep = sum(
				state.mass() for state in self.internal_states.itervalues())
			dry_mass_pre_timestep = all_submasses_pre_timestep.sum() - all_submasses_pre_timestep[
				self.water_index]

			# Calculate dry mass added in this timestep
			all_submass_diffs = sum(
				state.process_mass_diffs().sum(axis=0)
				for state in self.internal_states.itervalues())
			dry_mass_added = all_submass_diffs.sum() - all_submass_diffs[
				self.water_index]

			# Calculate dry mass after timestep
			dry_mass_post_timestep = dry_mass_pre_timestep + dry_mass_added

			# Log dry mass at t=0
			if not self.initial_mass_set:
				self.initial_mass_set = True
				self.initial_dry_mass = dry_mass_pre_timestep

			# End simulation once the mass of an average cell is added to the
			# current cell.
			current_media_id = self._external_states['Environment'].current_media_id
			if dry_mass_post_timestep - self.initial_dry_mass >= self.expectedDryMassIncreaseDict[
					current_media_id].asNumber(units.fg) * self.division_mass_multiplier:
				self._sim.cellCycleComplete()

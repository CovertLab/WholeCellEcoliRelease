"""
DnaSupercoiling listener

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/10/20
"""

from __future__ import division, absolute_import, print_function

import numpy as np

import wholecell.listeners.listener

class DnaSupercoiling(wholecell.listeners.listener.Listener):
	""" DnaSupercoiling """

	_name = 'DnaSupercoiling'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(DnaSupercoiling, self).__init__(*args, **kwargs)

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(DnaSupercoiling, self).initialize(sim, sim_data)

		self.uniqueMolecules = sim.internal_states['UniqueMolecules']
		self.relaxed_DNA_base_pairs_per_turn = sim_data.process.chromosome_structure.relaxed_DNA_base_pairs_per_turn

	# Allocate memory
	def allocate(self):
		super(DnaSupercoiling, self).allocate()

		self.segment_left_boundary_coordinates = np.array([], dtype=np.int64)
		self.segment_right_boundary_coordinates = np.array([], dtype=np.int64)
		self.segment_domain_indexes = np.array([], dtype=np.int32)
		self.segment_superhelical_densities = np.array([], dtype=np.float64)

	def update(self):
		# Get attributes of chromosomal segments
		chromosomal_segments = self.uniqueMolecules.container.objectsInCollection(
			'chromosomal_segment')
		boundary_coordinates, domain_indexes, linking_numbers = chromosomal_segments.attrs(
			'boundary_coordinates', 'domain_index', 'linking_number')

		# Get mask for segments with nonzero lengths
		segment_lengths = boundary_coordinates[:, 1] - boundary_coordinates[:, 0]

		assert np.all(segment_lengths >= 0)
		nonzero_length_mask = (segment_lengths > 0)

		self.segment_left_boundary_coordinates = boundary_coordinates[
			nonzero_length_mask, 0]
		self.segment_right_boundary_coordinates = boundary_coordinates[
			nonzero_length_mask, 1]
		self.segment_domain_indexes = domain_indexes[nonzero_length_mask]

		# Calculate superhelical densities
		linking_numbers_relaxed_DNA = (
			segment_lengths[nonzero_length_mask] / self.relaxed_DNA_base_pairs_per_turn)
		self.segment_superhelical_densities = np.divide(
			linking_numbers[nonzero_length_mask] - linking_numbers_relaxed_DNA,
			linking_numbers_relaxed_DNA)

	def tableCreate(self, tableWriter):
		tableWriter.set_variable_length_columns(
			'segment_left_boundary_coordinates',
			'segment_right_boundary_coordinates',
			'segment_domain_indexes',
			'segment_superhelical_densities',
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			segment_left_boundary_coordinates = self.segment_left_boundary_coordinates,
			segment_right_boundary_coordinates = self.segment_right_boundary_coordinates,
			segment_domain_indexes = self.segment_domain_indexes,
			segment_superhelical_densities = self.segment_superhelical_densities,
			)

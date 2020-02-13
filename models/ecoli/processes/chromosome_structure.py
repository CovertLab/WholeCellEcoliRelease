"""
ChromosomeStructure process

- Resolve physical collisions between molecules on the chromosome.
- TODO: Calculate the superhelical densities of each chromosomal segment

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/7/2020
"""

from __future__ import division, absolute_import, print_function

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences
from wholecell.listeners.listener import WriteMethod


class ChromosomeStructure(wholecell.processes.process.Process):
	""" ChromosomeStructure """

	_name = "ChromosomeStructure"

	# Constructor
	def __init__(self):
		super(ChromosomeStructure, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(ChromosomeStructure, self).initialize(sim, sim_data)

		# Load parameters
		self.RNA_sequences = sim_data.process.transcription.transcriptionSequences

		# Load bulk molecule views
		self.inactive_RNAPs = self.bulkMoleculeView(sim_data.moleculeIds.rnapFull)
		self.fragmentBases = self.bulkMoleculesView(
			[id_ + '[c]' for id_ in sim_data.moleculeGroups.fragmentNT_IDs])
		self.ppi = self.bulkMoleculeView(sim_data.moleculeIds.ppi)

		# Load unique molecule views
		self.active_replisomes = self.uniqueMoleculesView('active_replisome')
		self.active_RNAPs = self.uniqueMoleculesView('active_RNAP')
		self.RNAs = self.uniqueMoleculesView('RNA')


	def calculateRequest(self):
		if self.active_replisomes.total_counts()[0] > 0:
			# Request access to delete active RNAPs and RNAs
			self.active_RNAPs.request_access(self.EDIT_DELETE_ACCESS)
			self.RNAs.request_access(self.EDIT_DELETE_ACCESS)

	def evolveState(self):
		# If there are no active replisomes, return immediately
		if self.active_replisomes.total_counts()[0] == 0:
			return

		# Read attributes
		replisome_domain_indexes, replisome_coordinates = self.active_replisomes.attrs(
			'domain_index', 'coordinates')
		RNAP_domain_indexes, RNAP_coordinates, RNAP_directions, RNAP_unique_indexes = self.active_RNAPs.attrs(
			'domain_index', 'coordinates', 'direction', 'unique_index')
		RNA_TU_indexes, transcript_lengths, RNA_RNAP_indexes = self.RNAs.attrs(
			'TU_index', 'transcript_length', 'RNAP_index')

		# Build mask for RNAPs that should be removed
		RNAP_collision_mask = np.zeros_like(RNAP_domain_indexes, dtype=np.bool)

		# Loop through all domains with active replisomes
		for domain_index in np.unique(replisome_domain_indexes):
			domain_replisome_coordinates = replisome_coordinates[
				replisome_domain_indexes == domain_index]

			# Get mask for RNAPs on this domain that are "out of range"
			domain_collision_mask = np.logical_and.reduce((
				RNAP_domain_indexes == domain_index,
				RNAP_coordinates > domain_replisome_coordinates.min(),
				RNAP_coordinates < domain_replisome_coordinates.max()))

			RNAP_collision_mask[domain_collision_mask] = True

		# Build separate masks for head-on and co-directional collisions
		RNAP_headon_collision_mask = np.logical_and(
			RNAP_collision_mask,
			np.logical_xor(RNAP_directions, RNAP_coordinates > 0))
		RNAP_codirectional_collision_mask = np.logical_and(
			RNAP_collision_mask,
			np.logical_not(RNAP_headon_collision_mask))

		n_total_collisions = np.count_nonzero(RNAP_collision_mask)
		n_headon_collisions = np.count_nonzero(RNAP_headon_collision_mask)
		n_codirectional_collisions = np.count_nonzero(
			RNAP_codirectional_collision_mask)

		# Write values to listeners
		self.writeToListener(
			'RnapData', 'n_total_collisions', n_total_collisions)
		self.writeToListener(
			'RnapData', 'n_headon_collisions', n_headon_collisions)
		self.writeToListener(
			'RnapData', 'n_codirectional_collisions', n_codirectional_collisions)
		self.writeToListener(
			'RnapData', 'headon_collision_coordinates',
			RNAP_coordinates[RNAP_headon_collision_mask],
			writeMethod=WriteMethod.fill)
		self.writeToListener(
			'RnapData', 'codirectional_collision_coordinates',
			RNAP_coordinates[RNAP_codirectional_collision_mask],
			writeMethod=WriteMethod.fill)

		# Get mask for RNAs that are transcribed from removed RNAPs
		RNA_collision_mask = np.isin(
			RNA_RNAP_indexes, RNAP_unique_indexes[RNAP_collision_mask])

		# Remove RNAPs and RNAs that have collided with replisomes
		if n_total_collisions > 0:
			self.active_RNAPs.delByIndexes(np.where(RNAP_collision_mask)[0])
			self.RNAs.delByIndexes(np.where(RNA_collision_mask)[0])

			# Increment counts of inactive RNAPs
			self.inactive_RNAPs.countInc(n_total_collisions)

			# Get sequences of incomplete transcripts
			incomplete_sequence_lengths = transcript_lengths[
				RNA_collision_mask]
			n_initiated_sequences = np.count_nonzero(incomplete_sequence_lengths)

			if n_initiated_sequences > 0:
				incomplete_sequences = buildSequences(
					self.RNA_sequences,
					RNA_TU_indexes[RNA_collision_mask],
					np.zeros(n_total_collisions, dtype=np.int64),
					np.full(n_total_collisions, incomplete_sequence_lengths.max()))

				base_counts = np.zeros(4, dtype=np.int64)

				for sl, seq in zip(incomplete_sequence_lengths, incomplete_sequences):
					base_counts += np.bincount(seq[:sl], minlength=4)

				# Increment counts of fragment NTPs and phosphates
				self.fragmentBases.countsInc(base_counts)
				self.ppi.countInc(n_initiated_sequences)

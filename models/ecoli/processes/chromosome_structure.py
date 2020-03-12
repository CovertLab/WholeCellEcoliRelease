"""
ChromosomeStructure process

- Resolve collisions between molecules and replication forks on the chromosome.
- Remove and replicate promoters and motifs that are traversed by replisomes.
- Reset the boundaries and linking numbers of chromosomal segments.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/7/2020
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences
from six.moves import zip


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
		self.protein_sequences = sim_data.process.translation.translationSequences
		self.n_TUs = len(sim_data.process.transcription.rnaData)
		self.n_TFs = len(sim_data.process.transcription_regulation.tf_ids)
		self.n_amino_acids = len(sim_data.moleculeGroups.aaIDs)
		self.n_fragment_bases = len(sim_data.moleculeGroups.fragmentNT_IDs)
		replichore_lengths = sim_data.process.replication.replichore_lengths
		self.min_coordinates = -replichore_lengths[1]
		self.max_coordinates = replichore_lengths[0]
		self.relaxed_DNA_base_pairs_per_turn = sim_data.process.chromosome_structure.relaxed_DNA_base_pairs_per_turn
		self.terC_index = sim_data.process.chromosome_structure.terC_dummy_molecule_index
		
		# Load sim options
		self.calculate_superhelical_densities = sim._superhelical_density

		# Get placeholder value for chromosome domains without children
		self.no_child_place_holder = sim_data.process.replication.no_child_place_holder

		# Load bulk molecule views
		self.inactive_RNAPs = self.bulkMoleculeView(sim_data.moleculeIds.rnapFull)
		self.fragmentBases = self.bulkMoleculesView(
			[id_ + '[c]' for id_ in sim_data.moleculeGroups.fragmentNT_IDs])
		self.ppi = self.bulkMoleculeView(sim_data.moleculeIds.ppi)
		self.active_tfs = self.bulkMoleculesView(
			[x + "[c]" for x in sim_data.process.transcription_regulation.tf_ids])
		self.ribosome_30S_subunit = self.bulkMoleculeView(sim_data.moleculeIds.s30_fullComplex)
		self.ribosome_50S_subunit = self.bulkMoleculeView(sim_data.moleculeIds.s50_fullComplex)
		self.amino_acids = self.bulkMoleculesView(sim_data.moleculeGroups.aaIDs)
		self.water = self.bulkMoleculeView(sim_data.moleculeIds.water)

		# Load unique molecule views
		self.active_replisomes = self.uniqueMoleculesView('active_replisome')
		self.chromosome_domains = self.uniqueMoleculesView('chromosome_domain')
		self.active_RNAPs = self.uniqueMoleculesView('active_RNAP')
		self.RNAs = self.uniqueMoleculesView('RNA')
		self.active_ribosomes = self.uniqueMoleculesView('active_ribosome')
		self.oriCs = self.uniqueMoleculesView('oriC')
		self.full_chromosomes = self.uniqueMoleculesView('full_chromosome')
		self.promoters = self.uniqueMoleculesView('promoter')
		self.DnaA_boxes = self.uniqueMoleculesView('DnaA_box')
		
		if self.calculate_superhelical_densities:
			self.chromosomal_segments = self.uniqueMoleculesView(
				'chromosomal_segment')


	def calculateRequest(self):
		# Request access to delete active RNAPs, RNAs, ribosomes, chromosomal
		# segments, and DNA motifs
		self.active_RNAPs.request_access(self.EDIT_DELETE_ACCESS)
		self.RNAs.request_access(self.EDIT_DELETE_ACCESS)
		self.active_ribosomes.request_access(self.EDIT_DELETE_ACCESS)
		self.promoters.request_access(self.EDIT_DELETE_ACCESS)
		self.DnaA_boxes.request_access(self.EDIT_DELETE_ACCESS)
		
		if self.calculate_superhelical_densities:
			self.chromosomal_segments.request_access(self.EDIT_DELETE_ACCESS)

		# Request water to degrade polypeptides from removed ribosomes
		self.water.requestAll()


	def evolveState(self):
		# Read unique molecule attributes
		replisome_domain_indexes, replisome_coordinates, replisome_unique_indexes = self.active_replisomes.attrs(
			'domain_index', 'coordinates', 'unique_index')
		all_chromosome_domain_indexes, child_domains = self.chromosome_domains.attrs(
			'domain_index', 'child_domains')
		RNAP_domain_indexes, RNAP_coordinates, RNAP_directions, RNAP_unique_indexes = self.active_RNAPs.attrs(
			'domain_index', 'coordinates', 'direction', 'unique_index')

		origin_domain_indexes = self.oriCs.attr('domain_index')
		mother_domain_indexes = self.full_chromosomes.attr('domain_index')

		RNA_TU_indexes, transcript_lengths, RNA_RNAP_indexes, RNA_unique_indexes = self.RNAs.attrs(
			'TU_index', 'transcript_length', 'RNAP_index', 'unique_index')
		ribosome_protein_indexes, ribosome_peptide_lengths, ribosome_mRNA_indexes = self.active_ribosomes.attrs(
			'protein_index', 'peptide_length', 'mRNA_index')

		promoter_TU_indexes, promoter_domain_indexes, promoter_coordinates, promoter_bound_TFs = self.promoters.attrs(
			'TU_index', 'domain_index', 'coordinates', 'bound_TF')
		DnaA_box_domain_indexes, DnaA_box_coordinates, DnaA_box_bound = self.DnaA_boxes.attrs(
			'domain_index', 'coordinates', 'DnaA_bound')

		# Build dictionary of replisome coordinates with domain indexes as keys
		replisome_coordinates_from_domains = {
			domain_index: replisome_coordinates[replisome_domain_indexes == domain_index]
			for domain_index in np.unique(replisome_domain_indexes)
			}


		def get_removed_molecules_mask(domain_indexes, coordinates):
			"""
			Computes the boolean mask of unique molecules that should be
			removed based on the progression of the replication forks.
			"""
			mask = np.zeros_like(domain_indexes, dtype=np.bool)

			# Loop through all domains
			for domain_index in np.unique(domain_indexes):
				# Domain has active replisomes
				if domain_index in replisome_coordinates_from_domains:
					domain_replisome_coordinates = replisome_coordinates_from_domains[
						domain_index]

					# Get mask for molecules on this domain that are out of range
					domain_mask = np.logical_and.reduce((
						domain_indexes == domain_index,
						coordinates > domain_replisome_coordinates.min(),
						coordinates < domain_replisome_coordinates.max()))

				# Domain has no active replisomes
				else:
					# Domain has child domains (has finished replicating)
					if (child_domains[all_chromosome_domain_indexes == domain_index, 0]
							!= self.no_child_place_holder):
						# Remove all molecules on this domain
						domain_mask = (domain_indexes == domain_index)
					# Domain has not started replication
					else:
						continue

				mask[domain_mask] = True

			return mask


		# Build mask for molecules that should be removed
		removed_RNAPs_mask = get_removed_molecules_mask(
			RNAP_domain_indexes, RNAP_coordinates)
		removed_promoters_mask = get_removed_molecules_mask(
			promoter_domain_indexes, promoter_coordinates)
		removed_DnaA_boxes_mask = get_removed_molecules_mask(
			DnaA_box_domain_indexes, DnaA_box_coordinates)

		# Get attribute arrays of remaining RNAPs
		remaining_RNAPs_mask = np.logical_not(removed_RNAPs_mask)
		remaining_RNAP_domain_indexes = RNAP_domain_indexes[remaining_RNAPs_mask]
		remaining_RNAP_coordinates = RNAP_coordinates[remaining_RNAPs_mask]
		remaining_RNAP_unique_indexes = RNAP_unique_indexes[remaining_RNAPs_mask]

		# Build masks for head-on and co-directional collisions between RNAPs
		# and replication forks
		RNAP_headon_collision_mask = np.logical_and(
			removed_RNAPs_mask,
			np.logical_xor(RNAP_directions, RNAP_coordinates > 0))
		RNAP_codirectional_collision_mask = np.logical_and(
			removed_RNAPs_mask, np.logical_not(RNAP_headon_collision_mask))

		n_total_collisions = np.count_nonzero(removed_RNAPs_mask)
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
			RNAP_coordinates[RNAP_headon_collision_mask])
		self.writeToListener(
			'RnapData', 'codirectional_collision_coordinates',
			RNAP_coordinates[RNAP_codirectional_collision_mask])
		
		if self.calculate_superhelical_densities:
			# Get attributes of existing segments
			boundary_molecule_indexes, boundary_coordinates, segment_domain_indexes, linking_numbers = self.chromosomal_segments.attrs(
				'boundary_molecule_indexes', 'boundary_coordinates',
				'domain_index', 'linking_number')

			# Initialize new attributes of chromosomal segments
			all_new_boundary_molecule_indexes = np.empty((0, 2), dtype=np.int64)
			all_new_boundary_coordinates = np.empty((0, 2), dtype=np.int64)
			all_new_segment_domain_indexes = np.array([], dtype=np.int32)
			all_new_linking_numbers = np.array([], dtype=np.float64)
	
			for domain_index in np.unique(all_chromosome_domain_indexes):
				# Skip domains that have completed replication
				if np.all(domain_index < mother_domain_indexes):
					continue
	
				domain_spans_oriC = (domain_index in origin_domain_indexes)
				domain_spans_terC = (domain_index in mother_domain_indexes)

				# Get masks for segments and RNAPs in this domain
				segments_domain_mask = (segment_domain_indexes == domain_index)
				RNAP_domain_mask = (remaining_RNAP_domain_indexes == domain_index)

				# Parse attributes of segments in this domain
				boundary_molecule_indexes_this_domain = boundary_molecule_indexes[segments_domain_mask, :]
				boundary_coordinates_this_domain = boundary_coordinates[segments_domain_mask, :]
				linking_numbers_this_domain = linking_numbers[segments_domain_mask]

				# Parse attributes of remaining RNAPs in this domain
				new_molecule_coordinates_this_domain = remaining_RNAP_coordinates[RNAP_domain_mask]
				new_molecule_indexes_this_domain = remaining_RNAP_unique_indexes[RNAP_domain_mask]
	
				# Append coordinates and indexes of replisomes on this domain, if any
				if not domain_spans_oriC:
					replisome_domain_mask = (replisome_domain_indexes == domain_index)

					new_molecule_coordinates_this_domain = np.concatenate((
						new_molecule_coordinates_this_domain,
						replisome_coordinates[replisome_domain_mask]
						))
					new_molecule_indexes_this_domain = np.concatenate((
						new_molecule_indexes_this_domain,
						replisome_unique_indexes[replisome_domain_mask]
						))
	
				# Append coordinates and indexes of parent domain replisomes, if any
				if not domain_spans_terC:
					parent_domain_index = all_chromosome_domain_indexes[
						np.where(child_domains == domain_index)[0][0]]
					replisome_parent_domain_mask = (replisome_domain_indexes == parent_domain_index)

					new_molecule_coordinates_this_domain = np.concatenate((
						new_molecule_coordinates_this_domain,
						replisome_coordinates[replisome_parent_domain_mask]
						))
					new_molecule_indexes_this_domain = np.concatenate((
						new_molecule_indexes_this_domain,
						replisome_unique_indexes[replisome_parent_domain_mask]
						))
	
				# If there are no molecules left on this domain, continue
				if len(new_molecule_indexes_this_domain) == 0:
					continue

				# Calculate attributes of new segments
				new_segment_attrs = self._compute_new_segment_attributes(
					boundary_molecule_indexes_this_domain,
					boundary_coordinates_this_domain,
					linking_numbers_this_domain,
					new_molecule_indexes_this_domain,
					new_molecule_coordinates_this_domain,
					domain_spans_oriC, domain_spans_terC
					)
	
				# Append to existing array of new segment attributes
				all_new_boundary_molecule_indexes = np.vstack((
					all_new_boundary_molecule_indexes,
					new_segment_attrs["boundary_molecule_indexes"]))
				all_new_boundary_coordinates = np.vstack((
					all_new_boundary_coordinates,
					new_segment_attrs["boundary_coordinates"]))
				all_new_segment_domain_indexes = np.concatenate((
					all_new_segment_domain_indexes,
					np.full(len(new_segment_attrs["linking_numbers"]), domain_index,
						dtype=np.int32)))
				all_new_linking_numbers = np.concatenate((
					all_new_linking_numbers, new_segment_attrs["linking_numbers"]))
	
			# Delete all existing chromosomal segments
			self.chromosomal_segments.delByIndexes(
				np.arange(self.chromosomal_segments.total_counts()))
	
			# Add new chromosomal segments
			n_segments = len(all_new_linking_numbers)
			self.chromosomal_segments.moleculesNew(
				n_segments,
				boundary_molecule_indexes=all_new_boundary_molecule_indexes,
				boundary_coordinates=all_new_boundary_coordinates,
				domain_index=all_new_segment_domain_indexes,
				linking_number=all_new_linking_numbers,
				)


		# Get mask for RNAs that are transcribed from removed RNAPs
		removed_RNAs_mask = np.isin(
			RNA_RNAP_indexes, RNAP_unique_indexes[removed_RNAPs_mask])

		# Remove RNAPs and RNAs that have collided with replisomes
		if n_total_collisions > 0:
			self.active_RNAPs.delByIndexes(np.where(removed_RNAPs_mask)[0])
			self.RNAs.delByIndexes(np.where(removed_RNAs_mask)[0])

			# Increment counts of inactive RNAPs
			self.inactive_RNAPs.countInc(n_total_collisions)

			# Get sequences of incomplete transcripts
			incomplete_sequence_lengths = transcript_lengths[
				removed_RNAs_mask]
			n_initiated_sequences = np.count_nonzero(incomplete_sequence_lengths)

			if n_initiated_sequences > 0:
				incomplete_sequences = buildSequences(
					self.RNA_sequences,
					RNA_TU_indexes[removed_RNAs_mask],
					np.zeros(n_total_collisions, dtype=np.int64),
					np.full(n_total_collisions, incomplete_sequence_lengths.max()))

				base_counts = np.zeros(self.n_fragment_bases, dtype=np.int64)

				for sl, seq in zip(incomplete_sequence_lengths, incomplete_sequences):
					base_counts += np.bincount(seq[:sl], minlength=self.n_fragment_bases)

				# Increment counts of fragment NTPs and phosphates
				self.fragmentBases.countsInc(base_counts)
				self.ppi.countInc(n_initiated_sequences)

		# Get mask for ribosomes that are bound to removed mRNA molecules
		removed_ribosomes_mask = np.isin(
			ribosome_mRNA_indexes, RNA_unique_indexes[removed_RNAs_mask])
		n_removed_ribosomes = np.count_nonzero(removed_ribosomes_mask)

		# Remove ribosomes that are bound to removed mRNA molecules
		if n_removed_ribosomes > 0:
			self.active_ribosomes.delByIndexes(
				np.where(removed_ribosomes_mask)[0])

			# Increment counts of inactive ribosomal subunits
			self.ribosome_30S_subunit.countInc(n_removed_ribosomes)
			self.ribosome_50S_subunit.countInc(n_removed_ribosomes)

			# Get amino acid sequences of incomplete polypeptides
			incomplete_sequence_lengths = ribosome_peptide_lengths[
				removed_ribosomes_mask]
			n_initiated_sequences = np.count_nonzero(incomplete_sequence_lengths)

			if n_initiated_sequences > 0:
				incomplete_sequences = buildSequences(
					self.protein_sequences,
					ribosome_protein_indexes[removed_ribosomes_mask],
					np.zeros(n_removed_ribosomes, dtype=np.int64),
					np.full(n_removed_ribosomes, incomplete_sequence_lengths.max()))

				amino_acid_counts = np.zeros(
					self.n_amino_acids, dtype=np.int64)

				for sl, seq in zip(incomplete_sequence_lengths, incomplete_sequences):
					amino_acid_counts += np.bincount(
						seq[:sl], minlength=self.n_amino_acids)

				# Increment counts of free amino acids and decrease counts of
				# free water molecules
				self.amino_acids.countsInc(amino_acid_counts)
				self.water.countDec(
					incomplete_sequence_lengths.sum() - n_initiated_sequences)

		# Write to listener
		self.writeToListener(
			'RnapData', 'n_removed_ribosomes', n_removed_ribosomes)


		def get_replicated_motif_attributes(old_coordinates, old_domain_indexes):
			"""
			Computes the attributes of replicated motifs on the chromosome,
			given the old coordinates and domain indexes of the original motifs.
			"""
			# Coordinates are simply repeated
			new_coordinates = np.repeat(old_coordinates, 2)

			# Domain indexes are set to the child indexes of the original index
			new_domain_indexes = child_domains[
				np.array([np.where(all_chromosome_domain_indexes == idx)[0][0]
					for idx in old_domain_indexes]),
				:].flatten()

			return new_coordinates, new_domain_indexes


		# Replicate promoters
		n_new_promoters = 2*np.count_nonzero(removed_promoters_mask)

		if n_new_promoters > 0:
			# Delete original promoters
			self.promoters.delByIndexes(np.where(removed_promoters_mask)[0])

			# Add freed active tfs
			self.active_tfs.countsInc(
				promoter_bound_TFs[removed_promoters_mask, :].sum(axis=0))

			# Set up attributes for the replicated promoters
			promoter_TU_indexes_new = np.repeat(promoter_TU_indexes[removed_promoters_mask], 2)
			promoter_coordinates_new, promoter_domain_indexes_new = get_replicated_motif_attributes(
				promoter_coordinates[removed_promoters_mask],
				promoter_domain_indexes[removed_promoters_mask])

			# Add new promoters with new domain indexes
			self.promoters.moleculesNew(
				n_new_promoters,
				TU_index=promoter_TU_indexes_new,
				coordinates=promoter_coordinates_new,
				domain_index=promoter_domain_indexes_new,
				bound_TF=np.zeros((n_new_promoters, self.n_TFs), dtype=np.bool))


		# Replicate DnaA boxes
		n_new_DnaA_boxes = 2*np.count_nonzero(removed_DnaA_boxes_mask)

		if n_new_DnaA_boxes > 0:
			# Delete original DnaA boxes
			self.DnaA_boxes.delByIndexes(np.where(removed_DnaA_boxes_mask)[0])

			# Set up attributes for the replicated boxes
			DnaA_box_coordinates_new, DnaA_box_domain_indexes_new = get_replicated_motif_attributes(
				DnaA_box_coordinates[removed_DnaA_boxes_mask],
				DnaA_box_domain_indexes[removed_DnaA_boxes_mask])

			# Add new promoters with new domain indexes
			self.DnaA_boxes.moleculesNew(
				n_new_DnaA_boxes,
				coordinates=DnaA_box_coordinates_new,
				domain_index=DnaA_box_domain_indexes_new,
				DnaA_bound=np.zeros(n_new_DnaA_boxes, dtype=np.bool))


	def _compute_new_segment_attributes(self, old_boundary_molecule_indexes,
			old_boundary_coordinates, old_linking_numbers,
			new_molecule_indexes, new_molecule_coordinates,
			spans_oriC, spans_terC):
		# type: (np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, bool, bool) -> dict
		"""
		Calculates the updated attributes of chromosomal segments belonging to
		a specific chromosomal domain, given the previous and current
		coordinates of molecules bound to the chromosome.
		Args:
			old_boundary_molecule_indexes (np.ndarray, (N, 2)): Array of unique
				indexes of molecules that formed the boundaries of each
				chromosomal segment in the previous timestep.
			old_boundary_coordinates (np.ndarray, (N, 2)): Array of chromosomal
			 	coordinates of molecules that formed the boundaries of each
			 	chromosomal segment in the previous timestep.
			old_linking_numbers (np.ndarray, (N,)): Linking numbers of each
				chromosomal segment in the previous timestep.
			new_molecule_indexes (np.ndarray, (N,)): Unique indexes of all
				molecules bound to the domain at the current timestep.
			new_molecule_coordinates (np.ndarray, (N,)): Chromosomal
				coordinates of all molecules bound to the domain at the current
				timestep.
			spans_oriC (bool): True if the domain spans the origin.
			spans_terC (bool): True if the domain spans the terminus.
		Returns (wrapped as dict):
			boundary_molecule_indexes (np.ndarray, (M, 2)): Array of unique
				indexes of molecules that form the boundaries of new
				chromosomal segments.
			boundary_coordinates (np.ndarray, (M, 2)): Array of chromosomal
				coordinates of molecules that form the boundaries of new
				chromosomal segments.
			linking_numbers (np.ndarray, (M,)): Linking numbers of new
				chromosomal segments.
		"""
		# Sort old segment arrays by coordinates of left boundary
		old_coordinates_argsort = np.argsort(old_boundary_coordinates[:, 0])
		old_boundary_coordinates_sorted = old_boundary_coordinates[old_coordinates_argsort, :]
		old_boundary_molecule_indexes_sorted = old_boundary_molecule_indexes[old_coordinates_argsort, :]
		old_linking_numbers_sorted = old_linking_numbers[old_coordinates_argsort]

		# Sort new segment arrays by molecular coordinates
		new_coordinates_argsort = np.argsort(new_molecule_coordinates)
		new_molecule_coordinates_sorted = new_molecule_coordinates[new_coordinates_argsort]
		new_molecule_indexes_sorted = new_molecule_indexes[new_coordinates_argsort]

		# Domain does not span the origin
		if not spans_oriC:
			# A fragment spans oriC if two boundaries have opposite signs,
			# or both are equal to zero
			oriC_fragment_counts = np.count_nonzero(
				np.logical_not(np.logical_xor(
					old_boundary_coordinates_sorted[:, 0] < 0,
					old_boundary_coordinates_sorted[:, 1] > 0
					)))

			# if oriC fragment did not exist in the domain in the previous
			# timestep, add a dummy fragment that covers the origin with
			# linking number zero. This is done to generalize the
			# implementation of this method.
			if oriC_fragment_counts == 0:
				# Index of first segment where left boundary is nonnegative
				oriC_fragment_index = np.argmax(
					old_boundary_coordinates_sorted[:, 0] >= 0)

				# Get indexes of boundary molecules for this dummy segment
				oriC_fragment_boundary_molecule_indexes = np.array([
					old_boundary_molecule_indexes_sorted[oriC_fragment_index - 1, 1],
					old_boundary_molecule_indexes_sorted[oriC_fragment_index, 0]
					])

				# Insert dummy segment to array
				old_boundary_molecule_indexes_sorted = np.insert(
					old_boundary_molecule_indexes_sorted,
					oriC_fragment_index, oriC_fragment_boundary_molecule_indexes,
					axis=0)
				old_linking_numbers_sorted = np.insert(
					old_linking_numbers_sorted,
					oriC_fragment_index, 0)
			else:
				# There should not be more than one fragment that spans oriC
				assert oriC_fragment_counts == 1

		# Domain spans the terminus
		if spans_terC:
			# If the domain spans the terminus, dummy molecules are added to
			# each end of the chromosome s.t. the segment that spans terC is
			# split to two segments and we can maintain a linear representation
			# for the circular chromosome. These two segments are later
			# adjusted to have the same superhelical densities.
			new_molecule_coordinates_sorted = np.insert(
				new_molecule_coordinates_sorted,
				[0, len(new_molecule_coordinates_sorted)],
				[self.min_coordinates, self.max_coordinates]
				)

			new_molecule_indexes_sorted = np.insert(
				new_molecule_indexes_sorted,
				[0, len(new_molecule_indexes_sorted)], self.terC_index
				)

			# Add dummy molecule to old segments if they do not already exist
			if old_boundary_molecule_indexes_sorted[0, 0] != self.terC_index:
				old_boundary_molecule_indexes_sorted = np.vstack((
					np.array([self.terC_index, old_boundary_molecule_indexes_sorted[0, 0]]),
					old_boundary_molecule_indexes_sorted,
					np.array([old_boundary_molecule_indexes_sorted[-1, 1], self.terC_index])
					))
				old_linking_numbers_sorted = np.insert(
					old_linking_numbers_sorted,
					[0, len(old_linking_numbers_sorted)], 0
					)

		# Recalculate linking numbers of each segment after accounting for
		# boundary molecules that were removed in the current timestep
		linking_numbers_after_removal = []
		right_boundaries_retained = np.isin(
			old_boundary_molecule_indexes_sorted[:, 1],
			new_molecule_indexes_sorted)

		# Add up linking numbers of each segment until each retained boundary
		ln_this_fragment = 0.
		for retained, ln in zip(right_boundaries_retained, old_linking_numbers_sorted):
			ln_this_fragment += ln

			if retained:
				linking_numbers_after_removal.append(ln_this_fragment)
				ln_this_fragment = 0.

		# Number of segments should be equal to number of retained boundaries
		assert len(linking_numbers_after_removal) == right_boundaries_retained.sum()

		# Redistribute linking numbers of the two terC segments such that the
		# segments have same superhelical densities
		if spans_terC and np.count_nonzero(right_boundaries_retained) > 1:
			# Get molecule indexes of the boundaries of the two terC segments
			# left and right of terC
			retained_boundary_indexes = np.where(right_boundaries_retained)[0]
			left_segment_boundary_index = old_boundary_molecule_indexes_sorted[retained_boundary_indexes[0], 1]
			right_segment_boundary_index = old_boundary_molecule_indexes_sorted[retained_boundary_indexes[-2], 1]

			# Get mapping from molecule index to chromosomal coordinates
			molecule_index_to_coordinates = {
				index: coordinates for index, coordinates
				in zip(new_molecule_indexes_sorted, new_molecule_coordinates_sorted)
				}

			# Distribute linking number between two segments proportional to
			# the length of each segment
			left_segment_length = molecule_index_to_coordinates[left_segment_boundary_index] - self.min_coordinates
			right_segment_length = self.max_coordinates - molecule_index_to_coordinates[right_segment_boundary_index]
			full_segment_length = left_segment_length + right_segment_length
			full_linking_number = linking_numbers_after_removal[0] + linking_numbers_after_removal[-1]

			linking_numbers_after_removal[0] = full_linking_number * left_segment_length/full_segment_length
			linking_numbers_after_removal[-1] = full_linking_number * right_segment_length/full_segment_length

		# Get mask for molecules that already existed in the previous timestep
		existing_molecules_mask = np.isin(
			new_molecule_indexes_sorted, old_boundary_molecule_indexes_sorted)

		# Get numbers and lengths of new segments that each segment will be
		# split into
		segment_split_sizes = np.diff(np.where(existing_molecules_mask)[0])
		segment_lengths = np.diff(new_molecule_coordinates_sorted)

		assert len(segment_split_sizes) == len(linking_numbers_after_removal)

		# Calculate linking numbers of each segment after accounting for new
		# boundaries that were added
		new_linking_numbers = []
		i = 0

		for ln, size in zip(linking_numbers_after_removal, segment_split_sizes):
			if size == 1:
				new_linking_numbers.append(ln)
			else:
				# Split linking number proportional to length of segment
				total_length = segment_lengths[i:i + size].sum()
				new_linking_numbers.extend(
					list(ln*segment_lengths[i:i + size]/total_length)
					)
			i += size

		# Handle edge case where a domain was just initialized, and two
		# replisomes are bound to the origin
		if len(new_linking_numbers) == 0:
			new_linking_numbers = [0]

		# Build Mx2 array for boundary indexes and coordinates
		new_boundary_molecule_indexes = np.hstack((
			new_molecule_indexes_sorted[:-1, np.newaxis],
			new_molecule_indexes_sorted[1:, np.newaxis]
			))
		new_boundary_coordinates = np.hstack((
			new_molecule_coordinates_sorted[:-1, np.newaxis],
			new_molecule_coordinates_sorted[1:, np.newaxis]
			))
		new_linking_numbers = np.array(new_linking_numbers)

		# If domain does not span oriC, remove new segment that spans origin
		if not spans_oriC:
			oriC_fragment_mask = np.logical_not(np.logical_xor(
				new_boundary_coordinates[:, 0] < 0,
				new_boundary_coordinates[:, 1] > 0
				))

			assert oriC_fragment_mask.sum() == 1

			new_boundary_molecule_indexes = new_boundary_molecule_indexes[
				np.logical_not(oriC_fragment_mask), :]
			new_boundary_coordinates = new_boundary_coordinates[
				np.logical_not(oriC_fragment_mask), :]
			new_linking_numbers = new_linking_numbers[
				np.logical_not(oriC_fragment_mask)]

		return {
			"boundary_molecule_indexes": new_boundary_molecule_indexes,
			"boundary_coordinates": new_boundary_coordinates,
			"linking_numbers": new_linking_numbers
			}

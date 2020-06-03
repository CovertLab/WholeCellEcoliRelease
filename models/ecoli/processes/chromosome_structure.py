"""
ChromosomeStructure process

- Resolve collisions between molecules and replication forks on the chromosome.
- Remove and replicate promoters and motifs that are traversed by replisomes.
- TODO: Calculate the superhelical densities of each chromosomal segment

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/7/2020
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences


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
		self.promoters = self.uniqueMoleculesView('promoter')
		self.DnaA_boxes = self.uniqueMoleculesView('DnaA_box')


	def calculateRequest(self):
		# Request access to delete active RNAPs, RNAs, ribosomes, and DNA motifs
		self.active_RNAPs.request_access(self.EDIT_DELETE_ACCESS)
		self.RNAs.request_access(self.EDIT_DELETE_ACCESS)
		self.active_ribosomes.request_access(self.EDIT_DELETE_ACCESS)
		self.promoters.request_access(self.EDIT_DELETE_ACCESS)
		self.DnaA_boxes.request_access(self.EDIT_DELETE_ACCESS)

		# Request water to degrade polypeptides from removed ribosomes
		self.water.requestAll()

	def evolveState(self):
		# Read unique molecule attributes
		replisome_domain_indexes, replisome_coordinates = self.active_replisomes.attrs(
			'domain_index', 'coordinates')
		all_chromosome_domain_indexes, child_domains = self.chromosome_domains.attrs(
			'domain_index', 'child_domains')
		RNAP_domain_indexes, RNAP_coordinates, RNAP_directions, RNAP_unique_indexes = self.active_RNAPs.attrs(
			'domain_index', 'coordinates', 'direction', 'unique_index')
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

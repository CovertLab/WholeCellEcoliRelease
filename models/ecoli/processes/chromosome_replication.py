"""
Submodel for chromosome replication

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/12/2014
"""

from __future__ import division

import numpy as np
from itertools import izip

import wholecell.processes.process
from wholecell.utils.polymerize import (buildSequences, polymerize,
	computeMassIncrease)
from wholecell.utils import units


class ChromosomeReplication(wholecell.processes.process.Process):
	"""
	Performs initiation, elongation, and termination of active partial chromosomes
	that replicate the chromosome.
	"""

	_name = "ChromosomeReplication"

	# Constructor
	def __init__(self):
		super(ChromosomeReplication, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(ChromosomeReplication, self).initialize(sim, sim_data)

		# Load parameters
		self.criticalInitiationMass = sim_data.growthRateParameters.getDnaCriticalMass(
			sim_data.conditionToDoublingTime[sim_data.condition])
		self.getDnaCriticalMass = sim_data.growthRateParameters.getDnaCriticalMass
		self.nutrientToDoublingTime = sim_data.nutrientToDoublingTime
		self.replichore_lengths = sim_data.process.replication.replichore_lengths
		self.sequences = sim_data.process.replication.replication_sequences
		self.polymerized_dntp_weights = sim_data.process.replication.replicationMonomerWeights
		self.replication_coordinate = sim_data.process.transcription.rnaData[
			"replicationCoordinate"]
		self.D_period = sim_data.growthRateParameters.d_period.asNumber(
			units.s)

		# Create molecule views for replisome subunits, active replisomes,
		# origins of replication, chromosome domains, and free active TFs
		self.replisome_trimers = self.bulkMoleculesView(
			sim_data.moleculeGroups.replisome_trimer_subunits)
		self.replisome_monomers = self.bulkMoleculesView(
			sim_data.moleculeGroups.replisome_monomer_subunits)
		self.active_replisome = self.uniqueMoleculesView('active_replisome')
		self.oriCs = self.uniqueMoleculesView('originOfReplication')
		self.chromosome_domain = self.uniqueMoleculesView('chromosome_domain')
		self.active_tfs = self.bulkMoleculesView(
			[x + "[c]" for x in sim_data.process.transcription_regulation.tf_ids])

		# Create bulk molecule views for polymerization reaction
		self.dntps = self.bulkMoleculesView(sim_data.moleculeGroups.dNtpIds)
		self.ppi = self.bulkMoleculeView('PPI[c]')

		# Create molecules views for full chromosomes
		self.full_chromosome = self.uniqueMoleculesView("fullChromosome")

		# Create view for promoters and get total number of TF types
		self.promoters = self.uniqueMoleculesView("promoter")
		self.n_tf = len(sim_data.process.transcription_regulation.tf_ids)

		# Create view for DnaA boxes
		self.DnaA_boxes = self.uniqueMoleculesView("DnaA_box")

		# Get placeholder value for domains without children
		self.no_child_place_holder = sim_data.process.replication.no_child_place_holder

		self.base_elongation_rate = int(
			round(sim_data.growthRateParameters.dnaPolymeraseElongationRate.asNumber(
			units.nt / units.s)))
		self.replication = sim_data.process.replication

	def calculateRequest(self):
		# Get total count of existing oriC's
		n_oric = self.oriCs.total_counts()[0]

		# If there are no origins, return immediately
		if n_oric == 0:
			return

		# Get current cell mass
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)

		# Get critical initiation mass for current simulation environment
		current_media_id = self._external_states['Environment'].current_media_id
		self.criticalInitiationMass = self.getDnaCriticalMass(
			self.nutrientToDoublingTime[current_media_id])

		# Calculate mass per origin of replication, and compare to critical
		# initiation mass. If the cell mass has reached this critical mass,
		# the process will initiate a round of chromosome replication for each
		# origin of replication.
		massPerOrigin = cellMass / n_oric
		self.criticalMassPerOriC = massPerOrigin / self.criticalInitiationMass

		# If replication should be initiated, request subunits required for
		# building two replisomes per one origin of replication.
		if self.criticalMassPerOriC >= 1.0:
			self.replisome_trimers.requestIs(6*n_oric)
			self.replisome_monomers.requestIs(2*n_oric)

		# If there are no active forks return
		n_active_replisomes = self.active_replisome.total_counts()[0]
		if n_active_replisomes == 0:
			return

		# Get current locations of all replication forks
		active_replisomes = self.active_replisome.molecules_read_only()
		fork_coordinates = active_replisomes.attr("coordinates")
		sequence_length = np.abs(np.repeat(fork_coordinates, 2))

		rates = self.replication.make_elongation_rates(
			self.base_elongation_rate,
			self.timeStepSec())

		sequences = buildSequences(
			self.sequences,
			np.tile(np.arange(4), n_active_replisomes//2),
			sequence_length,
			rates)

		# Count number of each dNTP in sequences for the next timestep
		sequenceComposition = np.bincount(
			sequences[sequences != polymerize.PAD_VALUE], minlength=4)

		# If one dNTP is limiting then limit the request for the other three by
		# the same ratio
		dNtpsTotal = self.dntps.total_counts()
		maxFractionalReactionLimit = (np.fmin(1, dNtpsTotal / sequenceComposition)).min()

		# Request dNTPs
		self.dntps.requestIs(
			maxFractionalReactionLimit * sequenceComposition)

	def evolveState(self):
		## Module 1: Replication initiation
		# Get number of existing replisomes, oriCs, full chromosomes
		active_replisomes = self.active_replisome.molecules()
		n_active_replisomes = len(active_replisomes)

		oriCs = self.oriCs.molecules()
		n_oric = len(oriCs)

		full_chromosomes = self.full_chromosome.molecules()

		# Get existing chromosome domains and their attributes
		chromosome_domains = self.chromosome_domain.molecules()
		domain_index_existing_domain, child_domains = chromosome_domains.attrs(
			'domain_index', 'child_domains')

		# If there are no origins, return immediately
		if n_oric == 0:
			return

		# Get number of available replisome subunits
		n_replisome_trimers = self.replisome_trimers.counts()
		n_replisome_monomers = self.replisome_monomers.counts()

		# Initiate replication only when
		# 1) The cell has reached the critical mass per oriC
		# 2) There are enough replisome subunits to assemble two replisomes per
		# existing OriC.
		# Note that we assume asynchronous initiation does not happen.
		initiate_replication = (self.criticalMassPerOriC >= 1.0 and
			np.all(n_replisome_trimers == 6*n_oric) and
			np.all(n_replisome_monomers == 2*n_oric))

		# If all conditions are met, initiate a round of replication on every
		# origin of replication
		if initiate_replication:
			# Get attributes of existing oriCs and domains
			domain_index_existing_oric = oriCs.attr('domain_index')

			# Get indexes of the domains that would be getting child domains
			# (domains that contain an origin)
			new_parent_domains = np.where(np.in1d(domain_index_existing_domain,
				domain_index_existing_oric))[0]

			# Calculate counts of new replisomes and domains to add
			n_new_replisome = 2*n_oric
			n_new_domain = 2*n_oric

			# Calculate the domain indexes of new domains and oriC's
			max_domain_index = domain_index_existing_domain.max()
			domain_index_new = np.arange(
				max_domain_index + 1, max_domain_index + 2*n_oric + 1,
				dtype=np.int32)

			# Add new oriC's, and reset attributes of existing oriC's
			# All oriC's must be assigned new domain indexes
			oriCs.attrIs(domain_index=domain_index_new[:n_oric])
			self.oriCs.moleculesNew(
				n_oric, domain_index=domain_index_new[n_oric:])

			# Add and set attributes of newly created replisomes.
			# New replisomes inherit the domain indexes of the oriC's they
			# were initiated from. Two replisomes are formed per oriC, one on
			# the right replichore, and one on the left.
			coordinates_replisome = np.zeros(n_new_replisome, dtype=np.int64)
			right_replichore = np.tile(
				np.array([True, False], dtype=np.bool), n_oric)
			domain_index_new_replisome = np.repeat(
				domain_index_existing_oric, 2)

			self.active_replisome.moleculesNew(
				n_new_replisome,
				coordinates=coordinates_replisome,
				right_replichore=right_replichore,
				domain_index=domain_index_new_replisome)

			# Add and set attributes of new chromosome domains. All new domains
			# should have have no children domains.
			self.chromosome_domain.moleculesNew(
				n_new_domain,
				domain_index=domain_index_new,
				child_domains=np.full(
					(n_new_domain, 2), self.no_child_place_holder,
					dtype=np.int32))

			# Add new domains as children of existing domains
			child_domains[new_parent_domains] = domain_index_new.reshape(-1, 2)

			chromosome_domains.attrIs(
				child_domains=child_domains)

			# Decrement counts of replisome subunits
			self.replisome_trimers.countsDec(6*n_oric)
			self.replisome_monomers.countsDec(2*n_oric)

		# Write data from this module to a listener
		self.writeToListener("ReplicationData", "criticalMassPerOriC",
			self.criticalMassPerOriC)
		self.writeToListener("ReplicationData", "criticalInitiationMass",
			self.criticalInitiationMass.asNumber(units.fg))

		## Module 2: replication elongation
		# Get attributes of promoters
		promoters = self.promoters.molecules()

		TU_index, coordinates_promoters, domain_index_promoters, bound_TF = promoters.attrs(
			"TU_index", "coordinates", "domain_index", "bound_TF")

		# Write gene copy numbers to listener
		self.writeToListener(
			"RnaSynthProb", "gene_copy_number",
			np.bincount(TU_index, minlength=len(self.replication_coordinate)))

		# Get attributes of DnaA boxes
		DnaA_boxes = self.DnaA_boxes.molecules()

		coordinates_DnaA_boxes, domain_index_DnaA_boxes, DnaA_bound = DnaA_boxes.attrs(
			"coordinates", "domain_index", "DnaA_bound")

		# Write DnaA_box copy numbers to listener
		self.writeToListener(
			"ReplicationData", "total_DnaA_boxes", len(DnaA_boxes))
		self.writeToListener(
			"ReplicationData", "free_DnaA_boxes", np.logical_not(DnaA_bound).sum())

		# If no active replisomes are present, return immediately
		# Note: the new replication forks added in the previous module are not
		# elongated until the next timestep.
		if n_active_replisomes == 0:
			return

		# Get allocated counts of dNTPs
		dNtpCounts = self.dntps.counts()

		# Get attributes of existing replisomes
		domain_index_replisome, right_replichore, coordinates_replisome = active_replisomes.attrs(
			"domain_index", "right_replichore", "coordinates")
		sequence_length = np.abs(np.repeat(coordinates_replisome, 2))

		rates = self.replication.make_elongation_rates(
			self.base_elongation_rate,
			self.timeStepSec())

		# Build sequences to polymerize
		sequence_indexes = np.tile(np.arange(4), n_active_replisomes // 2)
		sequences = buildSequences(
			self.sequences,
			sequence_indexes,
			sequence_length,
			rates)

		# Use polymerize algorithm to quickly calculate the number of
		# elongations each fork catalyzes
		reactionLimit = dNtpCounts.sum()

		active_elongation_rates = self.elongation_rates[sequence_indexes]
		result = polymerize(
			sequences,
			dNtpCounts,
			reactionLimit,
			self.randomState,
			active_elongation_rates)

		sequenceElongations = result.sequenceElongation
		dNtpsUsed = result.monomerUsages

		# Compute mass increase for each elongated sequence
		mass_increase_dna = computeMassIncrease(
			sequences,
			sequenceElongations,
			self.polymerized_dntp_weights.asNumber(units.fg))

		# Compute masses that should be added to each replisome
		added_dna_mass = mass_increase_dna[0::2] + mass_increase_dna[1::2]

		# Update positions of each fork
		updated_length = sequence_length + sequenceElongations
		updated_coordinates = updated_length[0::2]

		# Reverse signs of fork coordinates on left replichore
		updated_coordinates[~right_replichore] = -updated_coordinates[~right_replichore]

		# Update attributes and submasses of replisomes
		active_replisomes.attrIs(coordinates = updated_coordinates)
		active_replisomes.add_submass_by_name("DNA", added_dna_mass)

		# Update counts of polymerized metabolites
		self.dntps.countsDec(dNtpsUsed)
		self.ppi.countInc(dNtpsUsed.sum())


		# Define function that identifies replicated DNA motifs
		def get_replicated_motif_mask(motif_coordinates, motif_domain_indexes):
			"""
			Computes a mask array for DNA motifs that should be replicated in
			this timestep, based on the old and new positions of replisomes.

			Args:
				motif_coordinates (ndarray): Replication coordinates of all
				existing motifs
				motif_domain_indexes (ndarray): Domain indexes of chromosome
				domains that each motif belongs to

			Returns: Mask array of motifs that should be replicated in this
			timestep
			"""
			# Initialize mask array
			replicated_motifs = np.zeros_like(motif_coordinates, dtype=np.bool)

			# Loop through all replisomes
			for (domain_index, rr, old_coord, new_coord) in izip(
					domain_index_replisome, right_replichore,
					coordinates_replisome, updated_coordinates):
				# Fork on right replichore
				if rr:
					coordinates_mask = np.logical_and(
						motif_coordinates >= old_coord,
						motif_coordinates < new_coord)

				# Fork on left replichore
				else:
					coordinates_mask = np.logical_and(
						motif_coordinates <= old_coord,
						motif_coordinates > new_coord)

				mask = np.logical_and(
					motif_domain_indexes == domain_index,
					coordinates_mask)

				replicated_motifs[mask] = True

			return replicated_motifs


		replicated_promoters = get_replicated_motif_mask(
			coordinates_promoters, domain_index_promoters)
		replicated_DnaA_boxes = get_replicated_motif_mask(
			coordinates_DnaA_boxes, domain_index_DnaA_boxes)

		# Get counts of replicated promoters and DnaA boxes
		n_new_promoters = 2*replicated_promoters.sum()
		n_new_DnaA_boxes = 2*replicated_DnaA_boxes.sum()

		# Handle replicated promoters
		if n_new_promoters > 0:
			# Delete original promoters
			promoters.delByIndexes(np.where(replicated_promoters)[0])

			# Add freed active tfs
			self.active_tfs.countsInc(
				bound_TF[replicated_promoters, :].sum(axis=0))

			# Set up attributes for the replicated promoters
			TU_index_new = np.repeat(TU_index[replicated_promoters], 2)
			coordinates_promoters_new = np.repeat(
				coordinates_promoters[replicated_promoters], 2)
			parent_domain_index_promoters = domain_index_promoters[replicated_promoters]

			domain_index_promoters_new = child_domains[
				np.array([np.where(domain_index_existing_domain == idx)[0][0]
					for idx in parent_domain_index_promoters]),
				:].flatten()

			# Add new promoters with new domain indexes
			self.promoters.moleculesNew(
				n_new_promoters,
				TU_index=TU_index_new,
				coordinates=coordinates_promoters_new,
				domain_index=domain_index_promoters_new,
				bound_TF=np.zeros((n_new_promoters, self.n_tf), dtype=np.bool))

		# Handle replicated DnaA boxes
		if n_new_DnaA_boxes > 0:
			# Delete original DnaA boxes
			DnaA_boxes.delByIndexes(np.where(replicated_DnaA_boxes)[0])

			# Set up attributes for the replicated DnaA boxes
			coordinates_DnaA_boxes_new = np.repeat(
				coordinates_DnaA_boxes[replicated_DnaA_boxes], 2)
			parent_domain_index_DnaA_boxes = domain_index_DnaA_boxes[
				replicated_DnaA_boxes]

			domain_index_DnaA_boxes_new = child_domains[
				np.array([np.where(domain_index_existing_domain == idx)[0][0]
					for idx in parent_domain_index_DnaA_boxes]),
				:].flatten()

			# Add new DnaA boxes with new domain indexes
			self.DnaA_boxes.moleculesNew(
				n_new_DnaA_boxes,
				coordinates=coordinates_DnaA_boxes_new,
				domain_index=domain_index_DnaA_boxes_new,
				DnaA_bound=np.zeros(n_new_DnaA_boxes, dtype=np.bool))

		## Module 3: replication termination
		# Determine if any forks have reached the end of their sequences. If
		# so, delete the replisomes and domains that were terminated.
		terminal_lengths = self.replichore_lengths[
			np.tile(np.arange(2), n_active_replisomes // 2)]
		terminated_replisomes = (np.abs(updated_coordinates) == terminal_lengths)

		# If any forks were terminated,
		if terminated_replisomes.sum() > 0:
			# Get domain indexes of terminated forks
			terminated_domains = np.unique(domain_index_replisome[terminated_replisomes])

			# Get attributes of existing domains and full chromosomes
			domain_index_domains, child_domains = chromosome_domains.attrs(
				"domain_index", "child_domains")
			domain_index_full_chroms = full_chromosomes.attr("domain_index")

			# Initialize array of replisomes and domains that should be deleted
			replisomes_to_delete = np.zeros_like(domain_index_replisome, dtype=np.bool)
			domains_to_delete = np.zeros_like(domain_index_domains, dtype=np.bool)

			# Count number of new full chromosomes that should be created
			n_new_chromosomes = 0

			# Initialize array for domain indexes of new full chromosomes
			domain_index_new_full_chroms = []

			for terminated_domain_index in terminated_domains:
				# Get all terminated replisomes in the terminated domain
				terminated_domain_matching_replisomes = np.logical_and(
					domain_index_replisome == terminated_domain_index,
					terminated_replisomes)

				# If both replisomes in the domain have terminated, we are
				# ready to split the chromosome and update the attributes.
				if terminated_domain_matching_replisomes.sum() == 2:
					# Tag replisomes and domains with the given domain index
					# for deletion
					replisomes_to_delete = np.logical_or(
						replisomes_to_delete,
						terminated_domain_matching_replisomes)

					domain_matching_domains = (
						domain_index_domains == terminated_domain_index)
					domains_to_delete = np.logical_or(
						domains_to_delete,
						domain_matching_domains)

					# Get child domains of deleted domain
					child_domains_this_domain = child_domains[
						np.where(domain_matching_domains)[0][0], :]

					# Modify domain index of one existing full chromosome to
					# index of first child domain
					domain_index_full_chroms[
						np.where(domain_index_full_chroms == terminated_domain_index)[0]
						] = child_domains_this_domain[0]

					# Increment count of new full chromosome
					n_new_chromosomes += 1

					# Append chromosome index of new full chromosome
					domain_index_new_full_chroms.append(child_domains_this_domain[1])

			# Delete terminated replisomes and domains
			active_replisomes.delByIndexes(np.where(replisomes_to_delete)[0])
			chromosome_domains.delByIndexes(np.where(domains_to_delete)[0])

			# Generate new full chromosome molecules
			if n_new_chromosomes > 0:
				self.full_chromosome.moleculesNew(
					n_new_chromosomes,
					division_time=[self.time() + self.D_period]*n_new_chromosomes,
					has_induced_division=[False] * n_new_chromosomes,
					domain_index=domain_index_new_full_chroms)

				# Reset domain index of existing chromosomes that have finished
				# replication
				full_chromosomes.attrIs(
					domain_index = domain_index_full_chroms)

			# Increment counts of replisome subunits
			self.replisome_trimers.countsInc(3*replisomes_to_delete.sum())
			self.replisome_monomers.countsInc(replisomes_to_delete.sum())


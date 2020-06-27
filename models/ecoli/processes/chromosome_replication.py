"""
Submodel for chromosome replication

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/12/2014
"""

from __future__ import absolute_import, division, print_function

import numpy as np

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

		self.max_time_step = sim_data.process.replication.max_time_step

		# Load parameters
		self.get_dna_critical_mass = sim_data.mass.get_dna_critical_mass
		self.criticalInitiationMass = self.get_dna_critical_mass(
			sim_data.conditionToDoublingTime[sim_data.condition])
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
		self.active_replisomes = self.uniqueMoleculesView('active_replisome')
		self.oriCs = self.uniqueMoleculesView('oriC')
		self.chromosome_domains = self.uniqueMoleculesView('chromosome_domain')

		# Create bulk molecule views for polymerization reaction
		self.dntps = self.bulkMoleculesView(sim_data.moleculeGroups.dNtpIds)
		self.ppi = self.bulkMoleculeView(sim_data.moleculeIds.ppi)

		# Create molecules views for full chromosomes
		self.full_chromosomes = self.uniqueMoleculesView('full_chromosome')

		# Get placeholder value for domains without children
		self.no_child_place_holder = sim_data.process.replication.no_child_place_holder

		self.basal_elongation_rate = int(
			round(sim_data.growthRateParameters.dnaPolymeraseElongationRate.asNumber(
			units.nt / units.s)))
		self.make_elongation_rates = sim_data.process.replication.make_elongation_rates


	def calculateRequest(self):
		# Get total count of existing oriC's
		n_oric = self.oriCs.total_count()

		# If there are no origins, return immediately
		if n_oric == 0:
			return

		# Get current cell mass
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)

		# Get critical initiation mass for current simulation environment
		current_media_id = self._external_states['Environment'].current_media_id
		self.criticalInitiationMass = self.get_dna_critical_mass(
			self.nutrientToDoublingTime[current_media_id])

		# Calculate mass per origin of replication, and compare to critical
		# initiation mass. If the cell mass has reached this critical mass,
		# the process will initiate a round of chromosome replication for each
		# origin of replication.
		massPerOrigin = cellMass / n_oric
		self.criticalMassPerOriC = massPerOrigin / self.criticalInitiationMass

		# If replication should be initiated, request subunits required for
		# building two replisomes per one origin of replication, and edit
		# access to oriC and chromosome domain attributes
		if self.criticalMassPerOriC >= 1.0:
			self.replisome_trimers.requestIs(6*n_oric)
			self.replisome_monomers.requestIs(2*n_oric)
			self.oriCs.request_access(self.EDIT_ACCESS)
			self.chromosome_domains.request_access(self.EDIT_ACCESS)

		# If there are no active forks return
		n_active_replisomes = self.active_replisomes.total_count()
		if n_active_replisomes == 0:
			return

		# Get current locations of all replication forks
		fork_coordinates = self.active_replisomes.attr("coordinates")
		sequence_length = np.abs(np.repeat(fork_coordinates, 2))

		self.elongation_rates = self.make_elongation_rates(
			self.randomState,
			len(self.sequences),
			self.basal_elongation_rate,
			self.timeStepSec())

		sequences = buildSequences(
			self.sequences,
			np.tile(np.arange(4), n_active_replisomes//2),
			sequence_length,
			self.elongation_rates)

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

		# Request access to relevant unique molecules
		self.full_chromosomes.request_access(self.EDIT_ACCESS)
		self.active_replisomes.request_access(self.EDIT_DELETE_ACCESS)

	def evolveState(self):
		## Module 1: Replication initiation
		# Get number of existing replisomes and oriCs
		n_active_replisomes = self.active_replisomes.total_count()
		n_oriC = self.oriCs.total_count()

		# If there are no origins, return immediately
		if n_oriC == 0:
			return

		# Get attributes of existing chromosome domains
		domain_index_existing_domain, child_domains = self.chromosome_domains.attrs(
			'domain_index', 'child_domains')

		# Get number of available replisome subunits
		n_replisome_trimers = self.replisome_trimers.counts()
		n_replisome_monomers = self.replisome_monomers.counts()

		# Initiate replication only when
		# 1) The cell has reached the critical mass per oriC
		# 2) There are enough replisome subunits to assemble two replisomes per
		# existing OriC.
		# Note that we assume asynchronous initiation does not happen.
		initiate_replication = (self.criticalMassPerOriC >= 1.0 and
			np.all(n_replisome_trimers == 6*n_oriC) and
			np.all(n_replisome_monomers == 2*n_oriC))

		# If all conditions are met, initiate a round of replication on every
		# origin of replication
		if initiate_replication:
			# Get attributes of existing oriCs and domains
			domain_index_existing_oric = self.oriCs.attr('domain_index')

			# Get indexes of the domains that would be getting child domains
			# (domains that contain an origin)
			new_parent_domains = np.where(np.in1d(domain_index_existing_domain,
				domain_index_existing_oric))[0]

			# Calculate counts of new replisomes and domains to add
			n_new_replisome = 2*n_oriC
			n_new_domain = 2*n_oriC

			# Calculate the domain indexes of new domains and oriC's
			max_domain_index = domain_index_existing_domain.max()
			domain_index_new = np.arange(
				max_domain_index + 1, max_domain_index + 2*n_oriC + 1,
				dtype=np.int32)

			# Add new oriC's, and reset attributes of existing oriC's
			# All oriC's must be assigned new domain indexes
			self.oriCs.attrIs(domain_index=domain_index_new[:n_oriC])
			self.oriCs.moleculesNew(
				n_oriC, domain_index=domain_index_new[n_oriC:])

			# Add and set attributes of newly created replisomes.
			# New replisomes inherit the domain indexes of the oriC's they
			# were initiated from. Two replisomes are formed per oriC, one on
			# the right replichore, and one on the left.
			coordinates_replisome = np.zeros(n_new_replisome, dtype=np.int64)
			right_replichore = np.tile(
				np.array([True, False], dtype=np.bool), n_oriC)
			domain_index_new_replisome = np.repeat(
				domain_index_existing_oric, 2)

			self.active_replisomes.moleculesNew(
				n_new_replisome,
				coordinates=coordinates_replisome,
				right_replichore=right_replichore,
				domain_index=domain_index_new_replisome)

			# Add and set attributes of new chromosome domains. All new domains
			# should have have no children domains.
			self.chromosome_domains.moleculesNew(
				n_new_domain,
				domain_index=domain_index_new,
				child_domains=np.full(
					(n_new_domain, 2), self.no_child_place_holder,
					dtype=np.int32))

			# Add new domains as children of existing domains
			child_domains[new_parent_domains] = domain_index_new.reshape(-1, 2)
			self.chromosome_domains.attrIs(child_domains=child_domains)

			# Decrement counts of replisome subunits
			self.replisome_trimers.countsDec(6*n_oriC)
			self.replisome_monomers.countsDec(2*n_oriC)

		# Write data from this module to a listener
		self.writeToListener("ReplicationData", "criticalMassPerOriC",
			self.criticalMassPerOriC)
		self.writeToListener("ReplicationData", "criticalInitiationMass",
			self.criticalInitiationMass.asNumber(units.fg))

		## Module 2: replication elongation
		# If no active replisomes are present, return immediately
		# Note: the new replication forks added in the previous module are not
		# elongated until the next timestep.
		if n_active_replisomes == 0:
			return

		# Get allocated counts of dNTPs
		dNtpCounts = self.dntps.counts()

		# Get attributes of existing replisomes
		domain_index_replisome, right_replichore, coordinates_replisome = self.active_replisomes.attrs(
			"domain_index", "right_replichore", "coordinates")
		sequence_length = np.abs(np.repeat(coordinates_replisome, 2))

		# Build sequences to polymerize
		sequence_indexes = np.tile(np.arange(4), n_active_replisomes // 2)

		sequences = buildSequences(
			self.sequences,
			sequence_indexes,
			sequence_length,
			self.elongation_rates)

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
		self.active_replisomes.attrIs(coordinates = updated_coordinates)
		self.active_replisomes.add_submass_by_name("DNA", added_dna_mass)

		# Update counts of polymerized metabolites
		self.dntps.countsDec(dNtpsUsed)
		self.ppi.countInc(dNtpsUsed.sum())


		## Module 3: replication termination
		# Determine if any forks have reached the end of their sequences. If
		# so, delete the replisomes and domains that were terminated.
		terminal_lengths = self.replichore_lengths[
			np.logical_not(right_replichore).astype(np.int64)]
		terminated_replisomes = (np.abs(updated_coordinates) == terminal_lengths)

		# If any forks were terminated,
		if terminated_replisomes.sum() > 0:
			# Get domain indexes of terminated forks
			terminated_domains = np.unique(domain_index_replisome[terminated_replisomes])

			# Get attributes of existing domains and full chromosomes
			domain_index_domains, child_domains = self.chromosome_domains.attrs(
				"domain_index", "child_domains")
			domain_index_full_chroms = self.full_chromosomes.attr("domain_index")

			# Initialize array of replisomes that should be deleted
			replisomes_to_delete = np.zeros_like(domain_index_replisome, dtype=np.bool)

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

					domain_mask = (
						domain_index_domains == terminated_domain_index)

					# Get child domains of deleted domain
					child_domains_this_domain = child_domains[
						np.where(domain_mask)[0][0], :]

					# Modify domain index of one existing full chromosome to
					# index of first child domain
					domain_index_full_chroms[
						np.where(domain_index_full_chroms == terminated_domain_index)[0]
						] = child_domains_this_domain[0]

					# Increment count of new full chromosome
					n_new_chromosomes += 1

					# Append chromosome index of new full chromosome
					domain_index_new_full_chroms.append(child_domains_this_domain[1])

			# Delete terminated replisomes
			self.active_replisomes.delByIndexes(np.where(replisomes_to_delete)[0])

			# Generate new full chromosome molecules
			if n_new_chromosomes > 0:
				self.full_chromosomes.moleculesNew(
					n_new_chromosomes,
					division_time=[self.time() + self.D_period]*n_new_chromosomes,
					has_triggered_division=[False] * n_new_chromosomes,
					domain_index=domain_index_new_full_chroms)

				# Reset domain index of existing chromosomes that have finished
				# replication
				self.full_chromosomes.attrIs(
					domain_index = domain_index_full_chroms)

			# Increment counts of replisome subunits
			self.replisome_trimers.countsInc(3*replisomes_to_delete.sum())
			self.replisome_monomers.countsInc(replisomes_to_delete.sum())

	def isTimeStepShortEnough(self, inputTimeStep, timeStepSafetyFraction):
		return inputTimeStep <= self.max_time_step

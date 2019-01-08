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
			sim_data.conditionToDoublingTime[sim_data.condition]
			)
		self.getDnaCriticalMass = sim_data.growthRateParameters.getDnaCriticalMass
		self.nutrientToDoublingTime = sim_data.nutrientToDoublingTime
		self.replichore_lengths = sim_data.process.replication.sequence_lengths[0::2]
		self.sequences = sim_data.process.replication.replication_sequences
		self.polymerized_dntp_weights = sim_data.process.replication.replicationMonomerWeights
		self.dnaPolyElngRate = int(
			round(sim_data.growthRateParameters.dnaPolymeraseElongationRate.asNumber(
			units.nt / units.s))
			)
		self.replication_coordinate = sim_data.process.transcription.rnaData[
			"replicationCoordinate"]
		self.D_period = sim_data.growthRateParameters.d_period.asNumber(
			units.s)

		# Create molecule views for replisome subunits, active replisomes,
		# origins of replication, and chromosome domains
		self.replisome_trimers = self.bulkMoleculesView(
			sim_data.moleculeGroups.replisome_trimer_subunits)
		self.replisome_monomers = self.bulkMoleculesView(
			sim_data.moleculeGroups.replisome_monomer_subunits)
		self.active_replisome = self.uniqueMoleculesView('active_replisome')
		self.oriCs = self.uniqueMoleculesView('originOfReplication')
		self.chromosome_domain = self.uniqueMoleculesView('chromosome_domain')

		# Create bulk molecule views for polymerization reaction
		self.dntps = self.bulkMoleculesView(sim_data.moleculeGroups.dNtpIds)
		self.ppi = self.bulkMoleculeView('PPI[c]')

		# Create bulk molecule view for gene copy number
		self.gene_copy_number = self.bulkMoleculesView(
			sim_data.process.transcription_regulation.geneCopyNumberColNames)

		# Create molecules views for full chromosomes
		self.full_chromosome = self.uniqueMoleculesView("fullChromosome")

		# Get placeholder value for domains without children
		self.no_child_place_holder = sim_data.process.replication.no_child_place_holder

	def calculateRequest(self):

		# Request all unique origins of replication
		self.oriCs.requestAll()

		# Get total count of existing oriC's
		n_oric = self.oriCs.total()[0]

		# Get current cell mass
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)

		# Get critical initiation mass for current simulation environment
		current_nutrients = self._external_states['Environment'].nutrients
		self.criticalInitiationMass = self.getDnaCriticalMass(
			self.nutrientToDoublingTime[current_nutrients])

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

		# Request all chromosome domains
		self.chromosome_domain.requestAll()

		# If there are no active forks return
		active_replisomes = self.active_replisome.allMolecules()
		n_active_replisomes = len(active_replisomes)
		if n_active_replisomes == 0:
			return

		# Request all replisomes and full chromosomes
		self.active_replisome.requestAll()
		self.full_chromosome.requestAll()

		# Get current locations of all replication forks
		fork_coordinates = active_replisomes.attr("coordinates")
		sequence_length = np.abs(np.repeat(fork_coordinates, 2))

		sequences = buildSequences(
			self.sequences,
			np.tile(np.arange(4), n_active_replisomes//2),
			sequence_length,
			self._dnaPolymeraseElongationRate()
			)

		# Count number of each dNTP in sequences for the next timestep
		sequenceComposition = np.bincount(
			sequences[sequences != polymerize.PAD_VALUE], minlength=4)

		# If one dNTP is limiting then limit the request for the other three by
		# the same ratio
		dNtpsTotal = self.dntps.total()
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

		# Get existing chromosome domains
		chromosome_domains = self.chromosome_domain.molecules()

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
			np.all(n_replisome_monomers == 2*n_oric)
			)

		# If all conditions are met, initiate a round of replication on every
		# origin of replication
		if initiate_replication:
			# Get attributes of existing oriCs and domains
			domain_index_existing_oric = oriCs.attr('domain_index')
			domain_index_existing_domain, child_domains = chromosome_domains.attrs(
				'domain_index', 'child_domains')

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
				dtype=np.int32
				)

			# Add new oriC's, replisomes, and domains
			oriCsNew = self.oriCs.moleculesNew(
				"originOfReplication",
				n_oric
				)
			active_replisomes_new = self.active_replisome.moleculesNew(
				"active_replisome",
				n_new_replisome
				)
			chromosome_domain_new = self.chromosome_domain.moleculesNew(
				"chromosome_domain",
				n_new_domain
				)

			# Set attributes of new oriC's, and reset attributes of existing
			# oriC's
            # All oriC's must be assigned new domain indexes
			oriCs.attrIs(
				domain_index=domain_index_new[:n_oric],
				)
			oriCsNew.attrIs(
				domain_index=domain_index_new[n_oric:],
				)

			# Calculate and set attributes of newly created replisomes.
			# New replisomes inherit the domain indexes of the oriC's they
			# were initiated from. Two replisomes are formed per oriC, one on
			# the right replichore, and one on the left.
			coordinates = np.zeros(n_new_replisome, dtype=np.int64)
			right_replichore = np.tile(
				np.array([True, False], dtype=np.bool), n_oric)
			domain_index_new_replisome = np.repeat(
				domain_index_existing_oric, 2)

			active_replisomes_new.attrIs(
				coordinates=coordinates,
				right_replichore=right_replichore,
				domain_index=domain_index_new_replisome,
				)

			# Set attributes of new chromosome domains. All new domains have
			# no children domains
			chromosome_domain_new.attrIs(
				domain_index=domain_index_new,
				child_domains=np.full(
					(n_new_domain, 2), self.no_child_place_holder,
					dtype=np.int32)
				)

			# Add new domains as children of existing domains
			child_domains[new_parent_domains] = domain_index_new.reshape(-1, 2)

			chromosome_domains.attrIs(
				child_domains=child_domains
				)

			# Decrement counts of replisome subunits
			self.replisome_trimers.countsDec(6*n_oric)
			self.replisome_monomers.countsDec(2*n_oric)

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
		domain_index_replisome, right_replichore, coordinates, massDiff_DNA = active_replisomes.attrs(
			"domain_index", "right_replichore", "coordinates", "massDiff_DNA"
			)
		sequence_length = np.abs(np.repeat(coordinates, 2))

		# Build sequences to polymerize
		sequences = buildSequences(
			self.sequences,
			np.tile(np.arange(4), n_active_replisomes // 2),
			sequence_length,
			self._dnaPolymeraseElongationRate()
			)

		# Use polymerize algorithm to quickly calculate the number of
		# elongations each fork catalyzes
		reactionLimit = dNtpCounts.sum()

		result = polymerize(
			sequences,
			dNtpCounts,
			reactionLimit,
			self.randomState
			)

		sequenceElongations = result.sequenceElongation
		dNtpsUsed = result.monomerUsages

		# Compute mass increase for each elongated sequence
		mass_increase_dna = computeMassIncrease(
			sequences,
			sequenceElongations,
			self.polymerized_dntp_weights.asNumber(units.fg)
			)

		# Compute masses that should be added to each replisome
		updatedMass = massDiff_DNA + mass_increase_dna[0::2] + mass_increase_dna[1::2]

		# Update positions of each fork
		updated_length = sequence_length + sequenceElongations
		updated_coordinates = updated_length[0::2]

		# Reverse signs of fork coordinates on left replichore
		updated_coordinates[~right_replichore] = -updated_coordinates[~right_replichore]

		active_replisomes.attrIs(
			coordinates = updated_coordinates,
			massDiff_DNA = updatedMass,
			)

		# Update counts of polymerized metabolites
		self.dntps.countsDec(dNtpsUsed)
		self.ppi.countInc(dNtpsUsed.sum())

		# Increment copy numbers of replicated genes
		new_gene_copies = np.zeros(len(self.replication_coordinate))

		for (rr, old_coord, new_coord) in izip(
				right_replichore, coordinates, updated_coordinates):
			# Fork on right replichore
			if rr:
				new_gene_copies[np.logical_and(
					self.replication_coordinate >= old_coord,
					self.replication_coordinate < new_coord
					)] += 1
			# Fork on left replichore
			else:
				new_gene_copies[np.logical_and(
					self.replication_coordinate <= old_coord,
					self.replication_coordinate > new_coord
					)] += 1

		self.gene_copy_number.countsInc(new_gene_copies)

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
				"domain_index", "child_domains"
				)
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
					terminated_replisomes
					)

				# If both replisomes in the domain have terminated, we are
				# ready to split the chromosome and update the attributes.
				if terminated_domain_matching_replisomes.sum() == 2:
					# Tag replisomes and domains with the given domain index
					# for deletion
					replisomes_to_delete = np.logical_or(
						replisomes_to_delete,
						terminated_domain_matching_replisomes
						)

					domain_matching_domains = (
							domain_index_domains == terminated_domain_index)
					domains_to_delete = np.logical_or(
						domains_to_delete,
						domain_matching_domains
						)

					# Get child domains of deleted domain
					child_domains_this_domain = child_domains[
						np.where(domain_matching_domains)[0], :]

					# Modify domain index of one existing full chromosome to
					# index of first child domain
					domain_index_full_chroms[
						np.where(domain_index_full_chroms == terminated_domain_index)[0]
						] = child_domains_this_domain[:, 0]

					# Increment count of new full chromosome
					n_new_chromosomes += 1

					# Append chromosome index of new full chromosome
					domain_index_new_full_chroms.append(child_domains_this_domain[:, 1])

			# Delete terminated replisomes and domains
			active_replisomes.delByIndexes(np.where(replisomes_to_delete)[0])
			chromosome_domains.delByIndexes(np.where(domains_to_delete)[0])

			# Generate new full chromosome molecules
			if n_new_chromosomes > 0:
				new_full_chromosome = self.full_chromosome.moleculesNew(
					"fullChromosome", n_new_chromosomes
					)
				new_full_chromosome.attrIs(
					division_time = [self.time() + self.D_period]*n_new_chromosomes,
					has_induced_division = [False]*n_new_chromosomes,
					domain_index = domain_index_new_full_chroms,
					)

				# Reset domain index of existing chromosomes that have finished
				# replication
				full_chromosomes.attrIs(
					domain_index = domain_index_full_chroms,
					)

			# Increment counts of replisome subunits
			self.replisome_trimers.countsInc(3*replisomes_to_delete.sum())
			self.replisome_monomers.countsInc(replisomes_to_delete.sum())


	def _dnaPolymeraseElongationRate(self):
		# Calculates elongation rate scaled by the time step
		return self.dnaPolyElngRate * self.timeStepSec()

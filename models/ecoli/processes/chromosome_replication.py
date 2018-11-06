"""
ChromosomeReplication

@author: Nick Ruggero
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
	""" ChromosomeReplication """

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
		self.sequenceLengths = sim_data.process.replication.sequence_lengths
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

		# Create unique molecule views for dna polymerases/replication forks
		# and origins of replication
		self.activeDnaPoly = self.uniqueMoleculesView('dnaPolymerase')
		self.oriCs = self.uniqueMoleculesView('originOfReplication')

		# Create bulk molecule views for polymerization reaction
		self.dntps = self.bulkMoleculesView(sim_data.moleculeGroups.dNtpIds)
		self.ppi = self.bulkMoleculeView('PPI[c]')

		# Create bulk molecule view for gene copy number
		self.gene_copy_number = self.bulkMoleculesView(
			sim_data.process.transcription_regulation.geneCopyNumberColNames)

		# Create molecules views for bulk and unique full chromosomes
		self.bulk_full_chromosome = self.bulkMoleculeView("CHROM_FULL[c]")
		self.unique_full_chromosome = self.uniqueMoleculesView("fullChromosome")

	def calculateRequest(self):

		# Request all unique origins of replication and replication forks
		self.oriCs.requestAll()

		# If there are no active forks return
		activeDnaPoly = self.activeDnaPoly.allMolecules()
		if len(activeDnaPoly) == 0:
			return

		# Request all active replication forks
		self.activeDnaPoly.requestAll()

		# Get sequences for all active forks
		sequenceIdx, sequenceLength = activeDnaPoly.attrs(
			'sequenceIdx', 'sequenceLength'
			)

		sequences = buildSequences(
			self.sequences,
			sequenceIdx,
			sequenceLength,
			self._dnaPolymeraseElongationRate()
			)

		# Count number of each dNTP in sequences for the next timestep
		sequenceComposition = np.bincount(
			sequences[sequences != polymerize.PAD_VALUE], minlength=4)

		# If one dNTP is limiting then limit the request for the other three by the same ratio
		dNtpsTotal = self.dntps.total()
		maxFractionalReactionLimit = (np.fmin(1, dNtpsTotal / sequenceComposition)).min()

		# Request dNTPs
		self.dntps.requestIs(
			maxFractionalReactionLimit * sequenceComposition)

	def evolveState(self):

		## Module 1: Replication initiation
		# Get number of active DNA polymerases and oriCs
		activeDnaPoly = self.activeDnaPoly.molecules()
		n_active_polymerase = len(activeDnaPoly)
		activePolymerasePresent = (n_active_polymerase > 0)
		oriCs = self.oriCs.molecules()
		n_oric = len(oriCs)
		n_chromosomes = self.bulk_full_chromosome.total()[0]

		# If there are no chromosomes and oriC's, return immediately
		if n_oric == 0 and n_chromosomes == 0:
			return

		# Get cell mass
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)

		# Get critical initiation mass for simulation medium environment
		current_nutrients = self._external_states['Environment'].nutrients
		self.criticalInitiationMass = self.getDnaCriticalMass(
			self.nutrientToDoublingTime[current_nutrients])

		# Calculate mass per origin of replication, and compare to critical
		# initiation mass. This is a rearrangement of the equation:
		# 	(Cell Mass)/(Number of origins) > Critical mass
		# If the above inequality holds true, initiate a round of chromosome
		# replication for every origin of replication
		massFactor = cellMass / self.criticalInitiationMass
		massPerOrigin = massFactor / n_oric
		replication_initiated = False

		# If conditions are true, initiate a round of replication on every
		# origin of replication
		if massPerOrigin >= 1.0:
			replication_initiated = True

			# Get replication round indexes of active DNA polymerases
			if activePolymerasePresent:
				replicationRound = activeDnaPoly.attr('replicationRound')
			else:
				# Set to -1 to set values for new polymerases to 0
				replicationRound = np.array([-1])

			# Get chromosome indexes of current oriCs
			chromosomeIndexOriC = oriCs.attr('chromosomeIndex')

			# Calculate number of new DNA polymerases required (4 per origin)
			n_new_polymerase = 4 * n_oric

			# Add new polymerases and oriC's
			activeDnaPolyNew = self.activeDnaPoly.moleculesNew(
				"dnaPolymerase",
				n_new_polymerase
				)
			oriCsNew = self.oriCs.moleculesNew(
				"originOfReplication",
				n_oric
				)

			# Calculate and set attributes of newly created polymerases
			sequenceIdx = np.tile(
				np.array([0, 1, 2, 3], dtype=np.int8),
				n_oric)
			sequenceLength = np.zeros(n_new_polymerase, dtype=np.int64)
			replicationRound = np.ones(n_new_polymerase, dtype=np.int64)*(
				replicationRound.max() + 1)

			# Polymerases inherit index of the OriC's they were initiated from
			chromosomeIndexPolymerase = np.repeat(chromosomeIndexOriC, 4)

			activeDnaPolyNew.attrIs(
				sequenceIdx=sequenceIdx,
				sequenceLength=sequenceLength,
				replicationRound=replicationRound,
				chromosomeIndex=chromosomeIndexPolymerase,
				)

			# Calculate and set attributes of newly created oriCs
            # New OriC's share the index of the old OriC's they were
            # replicated from
			oriCsNew.attrIs(chromosomeIndex=chromosomeIndexOriC)

		# Write data from this module to a listener
		self.writeToListener("ReplicationData", "criticalMassPerOriC",
			massPerOrigin)
		self.writeToListener("ReplicationData", "criticalInitiationMass",
			self.criticalInitiationMass.asNumber(units.fg))

		## Module 2: replication elongation
		# If no active polymerases are present, return immediately
		# Note: the new DNA polymerases activated in the previous module are
		# not elongated until the next timestep.
		if n_active_polymerase == 0:
			return

		# Build sequences to polymerize
		dNtpCounts = self.dntps.counts()
		sequenceIdx, sequenceLengths, massDiffDna = activeDnaPoly.attrs(
			'sequenceIdx', 'sequenceLength', 'massDiff_DNA'
			)

		sequences = buildSequences(
			self.sequences,
			sequenceIdx,
			sequenceLengths,
			self._dnaPolymeraseElongationRate()
			)

		# Use polymerize algorithm to quickly calculate the number of
		# elongations each "polymerase" catalyzes
		reactionLimit = dNtpCounts.sum()

		result = polymerize(
			sequences,
			dNtpCounts,
			reactionLimit,
			self.randomState
			)

		sequenceElongations = result.sequenceElongation
		dNtpsUsed = result.monomerUsages

		# Compute mass increase for each polymerizing chromosome
		massIncreaseDna = computeMassIncrease(
			sequences,
			sequenceElongations,
			self.polymerized_dntp_weights.asNumber(units.fg)
			)

		updatedMass = massDiffDna + massIncreaseDna

		# Update positions of each "polymerase"
		updatedLengths = sequenceLengths + sequenceElongations

		activeDnaPoly.attrIs(
			sequenceLength=updatedLengths,
			massDiff_DNA=updatedMass
			)

		# Update counts of polymerized metabolites
		self.dntps.countsDec(dNtpsUsed)
		self.ppi.countInc(dNtpsUsed.sum())

		# Increment copy numbers of replicated genes
		new_gene_copies = np.zeros(len(self.replication_coordinate))

		for (seq_idx, old_len, new_len) in izip(
				sequenceIdx, sequenceLengths, updatedLengths):
			# Fork on forward strand
			if seq_idx == 0:
				new_gene_copies[np.logical_and(
					self.replication_coordinate >= old_len,
					self.replication_coordinate < new_len
					)] += 1
			# Fork on reverse strand
			elif seq_idx == 1:
				new_gene_copies[np.logical_and(
					self.replication_coordinate <= -old_len,
					self.replication_coordinate > -new_len
					)] += 1

		self.gene_copy_number.countsInc(new_gene_copies)

		## Module 3: replication termination
		# Determine if any polymerases reached the end of their sequences. If
		# so, terminate replication and update the attributes of the remaining
		# polymerases and OriC's to reflect the new chromosome structure.
		terminalLengths = self.sequenceLengths[sequenceIdx]
		didTerminate = (updatedLengths == terminalLengths)

		# If any of the polymerases were terminated, check if all polymerases
		# initiated the same round as the terminated polymerases has already
		# been removed - if they have, update attributes of the remaining
		# polymerases and oriC's, and remove the polymerases.
		if didTerminate.sum() > 0:
			# Get attributes from active DNA polymerases and oriC's
			chromosomeIndexPolymerase, replicationRound = activeDnaPoly.attrs(
				'chromosomeIndex', 'replicationRound'
				)
			chromosomeIndexOriC = oriCs.attr('chromosomeIndex')

			# If new DNAPs were added in this timestep, append attributes of
			# these polymerases
			if replication_initiated:
				chromosomeIndexPolymeraseNew, replicationRoundNew = activeDnaPolyNew.attrs(
					'chromosomeIndex', 'replicationRound'
					)
				chromosomeIndexOriCNew = oriCsNew.attr('chromosomeIndex')

				chromosomeIndexPolymerase = np.append(
					chromosomeIndexPolymerase,
					chromosomeIndexPolymeraseNew)
				replicationRound = np.append(replicationRound,
					replicationRoundNew)
				chromosomeIndexOriC = np.append(chromosomeIndexOriC,
					chromosomeIndexOriCNew)

				didTerminate = np.pad(didTerminate, (0, n_new_polymerase),
					'constant')

			# Check that all terminated polymerases were initiated in the same
			# replication round
			assert np.unique(replicationRound[didTerminate]).size == 1

			# Get chromosome indexes of the terminated polymerases
			chromosomeIndexesTerminated = np.unique(chromosomeIndexPolymerase[didTerminate])
			newChromosomeIndex = chromosomeIndexPolymerase.max() + 1

			# Get replication round index of the terminated polymerases
			terminatedRound = replicationRound[didTerminate][0]

			# Count number of new full chromosomes that should be created
			n_new_chromosomes = 0

			# Keep track of polymerases that should be deleted
			polymerasesToDelete = np.zeros_like(chromosomeIndexPolymerase,
				dtype=np.bool)

			for chromosomeIndexTerminated in chromosomeIndexesTerminated:
				# Get all remaining active polymerases initiated in the same
				# replication round and in the given chromosome
				replicationRoundMatch = (replicationRound == terminatedRound)
				chromosomeMatch = (chromosomeIndexPolymerase == chromosomeIndexTerminated)
				remainingPolymerases = np.logical_and(
					replicationRoundMatch, chromosomeMatch)

				# Get all terminated polymerases in the given chromosome
				terminatedPolymerases = np.logical_and(
					didTerminate, chromosomeMatch)

				# If all active polymerases are terminated polymerases, we are
				# ready to split the chromosome and update the attributes.
				if remainingPolymerases.sum() == terminatedPolymerases.sum():

					# For each set of polymerases initiated in the same
					# replication round, update the chromosome indexes to a new
					# index for half of the polymerases.
					for roundIdx in np.arange(terminatedRound + 1, replicationRound.max() + 1):
						replicationRoundMatch = (replicationRound == roundIdx)
						polymerasesToSplit = np.logical_and(
							replicationRoundMatch, chromosomeMatch)

						n_matches = polymerasesToSplit.sum()

						# Number of polymerases initiated in a single round
						# must be a multiple of eight.
						assert n_matches % 8 == 0

						# Update the chromosome indexes for half of the polymerases
						secondHalfIdx = np.where(polymerasesToSplit)[0][(n_matches // 2):]
						chromosomeIndexPolymerase[secondHalfIdx] = newChromosomeIndex

					# Get oriC's in the chromosome getting divided
					chromosomeMatchOriC = (chromosomeIndexOriC == chromosomeIndexTerminated)
					n_matches = chromosomeMatchOriC.sum()

					# Number of OriC's in a dividing chromosome should be even
					assert n_matches % 2 == 0

					# Update the chromosome indexes for half of the OriC's
					secondHalfIdx = np.where(chromosomeMatchOriC)[0][(n_matches // 2):]
					chromosomeIndexOriC[secondHalfIdx] = newChromosomeIndex

					# Increment the new chromosome index in case another
					# chromosome needs to be split
					newChromosomeIndex += 1

					# Add terminated polymerases to the list to delete
					polymerasesToDelete = np.logical_or(polymerasesToDelete,
						terminatedPolymerases)

					# Increment count of new full chromosome
					n_new_chromosomes += 1

			# If new DNAPs were added in the same timestep, partition indexes
			# and reset attributes of old and new DNAPs separately
			if replication_initiated:
				# Reset chromosomeIndex for old DNAPs and oriC's
				activeDnaPoly.attrIs(
					chromosomeIndex=chromosomeIndexPolymerase[:n_active_polymerase])
				oriCs.attrIs(chromosomeIndex=chromosomeIndexOriC[:n_oric])

				# Reset chromosomeIndex for new DNAPs and oriC's
				activeDnaPolyNew.attrIs(
					chromosomeIndex=chromosomeIndexPolymerase[n_active_polymerase:])
				oriCsNew.attrIs(chromosomeIndex=chromosomeIndexOriC[n_oric:])

				# Delete terminated polymerases
				activeDnaPoly.delByIndexes(
					np.where(polymerasesToDelete[:n_active_polymerase])[0])

			else:
				# Reset chromosomeIndex for DNAPs and oriC's
				activeDnaPoly.attrIs(chromosomeIndex=chromosomeIndexPolymerase)
				oriCs.attrIs(chromosomeIndex=chromosomeIndexOriC)

				# Delete terminated polymerases
				activeDnaPoly.delByIndexes(np.where(polymerasesToDelete)[0])

			# Generate new full chromosome molecules
			if n_new_chromosomes > 0:
				self.bulk_full_chromosome.countInc(n_new_chromosomes)
				new_unique_full_chromosome = self.unique_full_chromosome.moleculesNew(
					"fullChromosome", n_new_chromosomes
					)
				new_unique_full_chromosome.attrIs(
					division_time = [self.time() + self.D_period]*n_new_chromosomes
					)


	def _dnaPolymeraseElongationRate(self):
		# Calculates elongation rate scaled by the time step
		return self.dnaPolyElngRate * self.timeStepSec()

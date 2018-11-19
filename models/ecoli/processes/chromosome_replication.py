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
	Performs initiation, elongation, and termination of active DNA polymerases
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

		# Create molecule views for replisome subunits, active replisomes, DNA
		# polymerases, and origins of replication
		self.replisome_trimers = self.bulkMoleculesView(
			sim_data.moleculeGroups.replisome_trimer_subunits)
		self.replisome_monomers = self.bulkMoleculesView(
			sim_data.moleculeGroups.replisome_monomer_subunits)
		self.activeReplisome = self.uniqueMoleculesView('activeReplisome')
		self.activeDnaPoly = self.uniqueMoleculesView('dnaPolymerase')
		self.oriCs = self.uniqueMoleculesView('originOfReplication')

		# Create bulk molecule views for polymerization reaction
		self.dntps = self.bulkMoleculesView(sim_data.moleculeGroups.dNtpIds)
		self.ppi = self.bulkMoleculeView('PPI[c]')

		# Create bulk molecule view for gene copy number
		self.gene_copy_number = self.bulkMoleculesView(
			sim_data.process.transcription_regulation.geneCopyNumberColNames)

		# Create molecules views for full chromosomes
		self.full_chromosome = self.uniqueMoleculesView("fullChromosome")


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

		# If there are no active forks return
		activeDnaPoly = self.activeDnaPoly.allMolecules()
		if len(activeDnaPoly) == 0:
			return

		# Request all active DNA polymerases and replisomes
		self.activeDnaPoly.requestAll()
		self.activeReplisome.requestAll()

		# Request all full chromosomes
		self.full_chromosome.requestAll()

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
		# Get number of active DNA polymerases, replisomes, and oriCs
		activeDnaPoly = self.activeDnaPoly.molecules()
		n_active_polymerase = len(activeDnaPoly)
		activePolymerasePresent = (n_active_polymerase > 0)

		activeReplisome = self.activeReplisome.molecules()
		n_active_replisome = len(activeReplisome)

		oriCs = self.oriCs.molecules()
		n_oric = len(oriCs)
		full_chromosomes = self.full_chromosome.molecules()
		n_chromosomes = len(full_chromosomes)

		# If there are no chromosomes and oriC's, return immediately
		if n_oric == 0 and n_chromosomes == 0:
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
			# Get replication round indexes of active DNA polymerases
			if activePolymerasePresent:
				replicationRound = activeDnaPoly.attr('replicationRound')
			else:
				# Set to -1 to set values for new polymerases to 0
				replicationRound = np.array([-1])

			# Get chromosome indexes of current oriCs
			chromosomeIndexOriC = oriCs.attr('chromosomeIndex')

			# Calculate number of new DNA polymerases and replisomes required
			n_new_polymerase = 4*n_oric
			n_new_replisome = 2*n_oric

			# Add new polymerases, oriC's, and replisomes
			activeDnaPolyNew = self.activeDnaPoly.moleculesNew(
				"dnaPolymerase",
				n_new_polymerase
				)
			oriCsNew = self.oriCs.moleculesNew(
				"originOfReplication",
				n_oric
				)
			activeReplisomeNew = self.activeReplisome.moleculesNew(
				"activeReplisome",
				n_new_replisome
			)

			# Calculate and set attributes of newly created polymerases
			# Polymerases inherit the chromosome indexes of the OriC's they
			# were initiated from
			sequenceIdx = np.tile(
				np.array([0, 1, 2, 3], dtype=np.int8), n_oric)
			sequenceLength = np.zeros(n_new_polymerase, dtype=np.int64)
			replicationRoundPolymerase = np.ones(n_new_polymerase,
				dtype=np.int64)*(replicationRound.max() + 1)
			chromosomeIndexPolymerase = np.repeat(chromosomeIndexOriC, 4)

			activeDnaPolyNew.attrIs(
				sequenceIdx=sequenceIdx,
				sequenceLength=sequenceLength,
				replicationRound=replicationRoundPolymerase,
				chromosomeIndex=chromosomeIndexPolymerase,
				)

			# Calculate and set attributes of newly created oriCs
            # New OriC's inherit the chromosome indexes of the old OriC's they
			# were replicated from.
			oriCsNew.attrIs(chromosomeIndex=chromosomeIndexOriC)

			# Calculate and set attributes of newly created replisomes.
			# New replisomes inherit the chromosome indexes of the OriC's they
			# were initiated from.
			replicationRoundReplisome = np.ones(n_new_replisome,
				dtype=np.int64)*(replicationRound.max() + 1)
			chromosomeIndexReplisome = np.repeat(chromosomeIndexOriC, 2)

			activeReplisomeNew.attrIs(
				replicationRound=replicationRoundReplisome,
				chromosomeIndex=chromosomeIndexReplisome,
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
		terminatedPolymerases = (updatedLengths == terminalLengths)

		# If any of the polymerases were terminated, check if all polymerases
		# initiated the same round as the terminated polymerases has already
		# been removed - if they have, update attributes of the remaining
		# polymerases and oriC's, and remove the polymerases.
		if terminatedPolymerases.sum() > 0:
			# Get attributes from active DNA polymerases and oriC's
			chromosomeIndexPolymerase, replicationRoundPolymerase = activeDnaPoly.attrs(
				'chromosomeIndex', 'replicationRound'
				)
			chromosomeIndexReplisome, replicationRoundReplisome = activeReplisome.attrs(
				'chromosomeIndex', 'replicationRound'
				)
			chromosomeIndexOriC = oriCs.attr('chromosomeIndex')
			chromosomeIndexFullChromosome = full_chromosomes.attr('chromosomeIndex')

			# If new DNAPs were added in this timestep, append attributes of
			# these polymerases
			if initiate_replication:
				chromosomeIndexPolymeraseNew, replicationRoundPolymeraseNew = activeDnaPolyNew.attrs(
					'chromosomeIndex', 'replicationRound'
					)
				chromosomeIndexReplisomeNew, replicationRoundReplisomeNew = activeReplisomeNew.attrs(
					'chromosomeIndex', 'replicationRound'
					)
				chromosomeIndexOriCNew = oriCsNew.attr('chromosomeIndex')

				chromosomeIndexPolymerase = np.append(
					chromosomeIndexPolymerase, chromosomeIndexPolymeraseNew
					)
				replicationRoundPolymerase = np.append(
					replicationRoundPolymerase, replicationRoundPolymeraseNew
					)
				chromosomeIndexReplisome = np.append(
					chromosomeIndexReplisome, chromosomeIndexReplisomeNew
					)
				replicationRoundReplisome = np.append(
					replicationRoundReplisome, replicationRoundReplisomeNew
					)
				chromosomeIndexOriC = np.append(
					chromosomeIndexOriC, chromosomeIndexOriCNew
					)

				terminatedPolymerases = np.pad(terminatedPolymerases,
					(0, n_new_polymerase), 'constant')

			# Check that all terminated polymerases were initiated in the same
			# replication round
			assert np.unique(replicationRoundPolymerase[terminatedPolymerases]).size == 1

			# Get chromosome indexes of the terminated polymerases
			chromosomeIndexesTerminated = np.unique(chromosomeIndexPolymerase[terminatedPolymerases])
			newChromosomeIndex = chromosomeIndexFullChromosome.max() + 1

			# Get replication round index of the terminated polymerases
			terminatedRound = replicationRoundPolymerase[terminatedPolymerases][0]

			# Initialize array of replisomes that need to be terminated
			terminatedReplisomes = np.zeros_like(chromosomeIndexReplisome, dtype=bool)

			# Count number of new full chromosomes that should be created
			n_new_chromosomes = 0

			# Initialize chromosome indexes of new full chromosomes
			chromosomeIndexFullChromosomeNew = []

			# Keep track of polymerases that should be deleted
			polymerasesToDelete = np.zeros_like(chromosomeIndexPolymerase,
				dtype=np.bool)

			for chromosomeIndexTerminated in chromosomeIndexesTerminated:
				# Get all remaining active polymerases initiated in the same
				# replication round and in the given chromosome
				replicationRoundMatchPolymerase = (
						replicationRoundPolymerase == terminatedRound)
				chromosomeMatchPolymerase = (
						chromosomeIndexPolymerase == chromosomeIndexTerminated)
				remainingPolymerasesChromosome = np.logical_and(
					replicationRoundMatchPolymerase, chromosomeMatchPolymerase)

				# Get all terminated polymerases in the given chromosome
				terminatedPolymerasesChromosome = np.logical_and(
					terminatedPolymerases, chromosomeMatchPolymerase)

				# Get all active replisomes in the given chromosome
				chromosomeMatchReplisome = (
						chromosomeIndexReplisome == chromosomeIndexTerminated)

				# If all active polymerases are terminated polymerases, we are
				# ready to split the chromosome and update the attributes.
				if remainingPolymerasesChromosome.sum() == terminatedPolymerasesChromosome.sum():

					# For each set of polymerases/replisomes initiated in the
					# same replication round, update the chromosome indexes to
					# a new index for half of the polymerases/replisomes.
					for roundIdx in np.arange(terminatedRound + 1,
							replicationRoundPolymerase.max() + 1):

						replicationRoundMatchPolymerase = (
								replicationRoundPolymerase == roundIdx)
						polymerasesToSplit = np.logical_and(
							replicationRoundMatchPolymerase,
							chromosomeMatchPolymerase)
						n_matches_polymerase = polymerasesToSplit.sum()

						replicationRoundMatchReplisome = (
								replicationRoundReplisome == roundIdx)
						replisomesToSplit = np.logical_and(
							replicationRoundMatchReplisome,
							chromosomeMatchReplisome)
						n_matches_replisome = replisomesToSplit.sum()

						# Number of polymerases/replisomes initiated in a
						# single round must be a multiple of eight/four.
						assert n_matches_polymerase % 8 == 0
						assert n_matches_replisome % 4 == 0

						# Update the chromosome indexes for half of the polymerases
						secondHalfIdxPolymerase = np.where(
							polymerasesToSplit)[0][(n_matches_polymerase // 2):]
						chromosomeIndexPolymerase[secondHalfIdxPolymerase] = newChromosomeIndex

						secondHalfIdxReplisome = np.where(
							replisomesToSplit)[0][(n_matches_replisome // 2):]
						chromosomeIndexReplisome[secondHalfIdxReplisome] = newChromosomeIndex

					# Get oriC's in the chromosome getting divided
					chromosomeMatchOriC = (chromosomeIndexOriC == chromosomeIndexTerminated)
					n_matches = chromosomeMatchOriC.sum()

					# Number of OriC's in a dividing chromosome should be even
					assert n_matches % 2 == 0

					# Update the chromosome indexes for half of the OriC's
					secondHalfIdx = np.where(chromosomeMatchOriC)[0][(n_matches // 2):]
					chromosomeIndexOriC[secondHalfIdx] = newChromosomeIndex

					# Add replisomes with the same replication round and
					# chromosome index as the terminating DNAPs to the list
					# of replisomes to terminate.
					replicationRoundMatchReplisome = (
						replicationRoundReplisome == terminatedRound)

					terminatedReplisomes = np.logical_or(terminatedReplisomes,
						np.logical_and(replicationRoundMatchReplisome,
						chromosomeMatchReplisome)
						)

					# Add terminated polymerases to the list to delete
					polymerasesToDelete = np.logical_or(polymerasesToDelete,
						terminatedPolymerases)

					# Increment count of new full chromosome
					n_new_chromosomes += 1

					# Append chromosome index of new full chromosome
					chromosomeIndexFullChromosomeNew.append(newChromosomeIndex)

					# Increment the new chromosome index in case another
					# chromosome needs to be split
					newChromosomeIndex += 1

			# If new DNAPs and replisomes were added in the same timestep,
			# partition indexes and reset attributes of old and new DNAPs and
			# replisomes separately
			if initiate_replication:
				# Reset chromosomeIndex for old DNAPs, replisomes and oriC's
				activeDnaPoly.attrIs(
					chromosomeIndex=chromosomeIndexPolymerase[:n_active_polymerase]
					)
				activeReplisome.attrIs(
					chromosomeIndex=chromosomeIndexReplisome[:n_active_replisome]
					)
				oriCs.attrIs(chromosomeIndex=chromosomeIndexOriC[:n_oric])

				# Reset chromosomeIndex for new DNAPs, replisomes and oriC's
				activeDnaPolyNew.attrIs(
					chromosomeIndex=chromosomeIndexPolymerase[n_active_polymerase:]
					)
				activeReplisomeNew.attrIs(
					chromosomeIndex=chromosomeIndexReplisome[n_active_replisome]
					)
				oriCsNew.attrIs(chromosomeIndex=chromosomeIndexOriC[n_oric:])

				# Delete terminated polymerases
				activeDnaPoly.delByIndexes(
					np.where(polymerasesToDelete[:n_active_polymerase])[0]
					)
				activeReplisome.delByIndexes(
					np.where(terminatedReplisomes[:n_active_replisome])[0]
					)

			else:
				# Reset chromosomeIndex for DNAPs, replisomes, and oriC's
				activeDnaPoly.attrIs(chromosomeIndex=chromosomeIndexPolymerase)
				activeReplisome.attrIs(chromosomeIndex=chromosomeIndexReplisome)
				oriCs.attrIs(chromosomeIndex=chromosomeIndexOriC)

				# Delete terminated polymerases and replisomes
				activeDnaPoly.delByIndexes(np.where(polymerasesToDelete)[0])
				activeReplisome.delByIndexes(np.where(terminatedReplisomes)[0])

			# Generate new full chromosome molecules
			if n_new_chromosomes > 0:
				new_full_chromosome = self.full_chromosome.moleculesNew(
					"fullChromosome", n_new_chromosomes
					)
				new_full_chromosome.attrIs(
					division_time = [self.time() + self.D_period]*n_new_chromosomes,
					chromosomeIndex = chromosomeIndexFullChromosomeNew,
					)

			# Increment counts of replisome subunits
			self.replisome_trimers.countsInc(3*terminatedReplisomes.sum())
			self.replisome_monomers.countsInc(terminatedReplisomes.sum())


	def _dnaPolymeraseElongationRate(self):
		# Calculates elongation rate scaled by the time step
		return self.dnaPolyElngRate * self.timeStepSec()

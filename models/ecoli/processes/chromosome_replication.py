"""
ChromosomeReplication

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/12/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.polymerize import buildSequences, polymerize, computeMassIncrease
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
		self.dnaPolyElngRate = int(round(sim_data.growthRateParameters.dnaPolymeraseElongationRate.asNumber(units.nt / units.s)))

		# Create unique molecule views for dna polymerases/replication forks and origins of replication
		self.activeDnaPoly = self.uniqueMoleculesView('dnaPolymerase')
		self.oriCs = self.uniqueMoleculesView('originOfReplication')

		# Create bulk molecule views for polymerization reaction
		self.dntps = self.bulkMoleculesView(sim_data.moleculeGroups.dNtpIds)
		self.ppi = self.bulkMoleculeView('PPI[c]')
		self.partialChromosomes = self.bulkMoleculesView(sim_data.moleculeGroups.partialChromosome)

		# Create bulk molecules view for full chromosome
		self.full_chromosome = self.bulkMoleculeView("CHROM_FULL[c]")

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
		sequenceComposition = np.bincount(sequences[sequences != polymerize.PAD_VALUE], minlength = 4)

		# If one dNTP is limiting then limit the request for the other three by the same ratio
		dNtpsTotal = self.dntps.total()
		maxFractionalReactionLimit = (np.fmin(1, dNtpsTotal/sequenceComposition)).min()

		# Request dNTPs
		self.dntps.requestIs(
			maxFractionalReactionLimit * sequenceComposition
			)

	def evolveState(self):

		## Module 1: Replication initiation
		# Get number of active DNA polymerases and oriCs
		activeDnaPoly = self.activeDnaPoly.molecules()
		activePolymerasePresent = (len(activeDnaPoly) > 0)
		oriCs = self.oriCs.molecules()
		n_oric = len(oriCs)
		n_chromosomes = self.full_chromosome.total()[0]

		# If there are no chromosomes and oriC's, return immediately
		if n_oric == 0 and n_chromosomes == 0:
			return

		# Get cell mass
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)

		# Get critical initiation mass for simulation medium environment
		current_nutrients = self._external_states['Environment'].nutrients
		self.criticalInitiationMass = self.getDnaCriticalMass(self.nutrientToDoublingTime[current_nutrients])

		# Calculate mass per origin of replication, and compare to critical
		# initiation mass. This is a rearrangement of the equation:
		# 	(Cell Mass)/(Number of origins) > Critical mass
		# If the above inequality holds true, initiate a round of chromosome
		# replication for every origin of replication
		massFactor = cellMass / self.criticalInitiationMass
		massPerOrigin = massFactor / n_oric

		# If conditions are true, initiate a round of replication on every
		# origin of replication
		if massPerOrigin >= 1.0 and self.partialChromosomes.counts().sum() == 0:

			# Get replication round indexes of active DNA polymerases
			if activePolymerasePresent:
				replicationRound = activeDnaPoly.attr('replicationRound')
			else:
				replicationRound = np.array([-1])  # Set to -1 to set values for new polymerases to 0

			# Calculate number of new DNA polymerases required per origin
			n_new_polymerase = 4*n_oric

			# Initialize 4 new polymerases per origin
			activeDnaPoly_new = self.activeDnaPoly.moleculesNew(
				"dnaPolymerase",
				n_new_polymerase
				)

			# Generate a new origin for each old origin
			self.oriCs.moleculesNew("originOfReplication", n_oric)

			# Calculate and set attributes of newly created polymerases
			sequenceIdx = np.tile(np.array([0, 1, 2, 3], dtype=np.int8), n_oric)
			sequenceLength = np.zeros(n_new_polymerase, dtype=np.int8)
			replicationRound = np.ones(n_new_polymerase, dtype=np.int8)*(replicationRound.max() + 1)
			chromosomeIndex = np.zeros(n_new_polymerase, dtype=np.int8)
			# TODO: check if this is right
			chromosomeIndex[(n_new_polymerase//n_chromosomes):] = 1

			activeDnaPoly_new.attrIs(
				sequenceIdx = sequenceIdx,
				sequenceLength = sequenceLength,
				replicationRound = replicationRound,
				chromosomeIndex = chromosomeIndex,
				)

		# Write data from this module to a listener
		self.writeToListener("ReplicationData", "criticalMassPerOriC", massPerOrigin)
		self.writeToListener("ReplicationData", "criticalInitiationMass", self.criticalInitiationMass.asNumber(units.fg))

		## Module 2: replication elongation
		# If no active polymerases are present, return immediately
		# Note: the new DNA polymerases activated in the previous module are
		# not elongated until the next timestep.
		if len(activeDnaPoly) == 0:
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
			sequenceLength = updatedLengths,
			massDiff_DNA = updatedMass
			)

		# Update counts of polymerized metabolites
		self.dntps.countsDec(dNtpsUsed)
		self.ppi.countInc(dNtpsUsed.sum())

		## Module 3: replication termination
		# Determine if any polymerases reached the end of their sequences. If
		# so, terminate replication and update the attributes of the remaining
		# polymerases to reflect the new chromosome structure
		terminalLengths = self.sequenceLengths[sequenceIdx]
		didTerminate = (updatedLengths == terminalLengths)

		terminatedChromosomes = np.bincount(
			sequenceIdx[didTerminate],
			minlength = self.sequences.shape[0]
			)

		# If any forks terminated, update the chromosome index of the next
		# oldest replication fork
		if didTerminate.sum() > 0:
			sequenceIdx, sequenceLength, chromosomeIndex, replicationRound = activeDnaPoly.attrs(
				'sequenceIdx', 'sequenceLength', 'chromosomeIndex', 'replicationRound'
				)

			# Search for next oldest replication fork
			# Criteria: The fork must have the same sequence index as the
			# terminating fork, a replication round index that is one greater
			# than the terminating fork, and the same chromosome index as the
			# terminating fork (so you don't flip the same upstream fork
			# twice).
			sequenceIdxMatch = np.zeros(didTerminate.shape)
			replicationRoundMatch = np.zeros(didTerminate.shape)
			chromosomeIndexMatch = np.zeros(didTerminate.shape)

			for idx in np.where(didTerminate)[0]:
				sequenceIdxMatch = np.logical_or(sequenceIdxMatch, sequenceIdx[idx] == sequenceIdx)
				replicationRoundMatch = np.logical_or(replicationRoundMatch, replicationRound[idx] + 1 == replicationRound)
				chromosomeIndexMatch = np.logical_or(chromosomeIndexMatch, chromosomeIndex[idx] == chromosomeIndex)

			# Potential matches for chromosome index to fli
			potentialMatches = np.logical_and.reduce((sequenceIdxMatch, replicationRoundMatch, chromosomeIndexMatch))

			# Take the first two. If there are more they are degenerate.
			toIndex = np.where(potentialMatches)[0][:2]

			# Changes 0 -> 1 and 1 -> 0
			toFlip = chromosomeIndex[toIndex]
			chromosomeIndex[toIndex] = np.abs(toFlip - 1)

			# Resets chromosomeIndex
			activeDnaPoly.attrIs(chromosomeIndex = chromosomeIndex)

		# Delete terminated "polymerases"
		activeDnaPoly.delByIndexes(np.where(didTerminate)[0])

		# Update counts of newly created chromosome halves. These will be
		# "stitched" together in the ChromosomeFormation process
		self.partialChromosomes.countsInc(terminatedChromosomes)


	def _dnaPolymeraseElongationRate(self):
		# Calculates elongation rate scaled by the time step
		return self.dnaPolyElngRate * self.timeStepSec()

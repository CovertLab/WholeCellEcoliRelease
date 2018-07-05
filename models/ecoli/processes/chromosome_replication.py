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
		self.chromosomeHalves = self.bulkMoleculesView(sim_data.moleculeGroups.partialChromosome)

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

	# Calculate temporal evolution
	def evolveState(self):
		# Set critical initiaion mass for simulation medium environment
		current_nutrients = self._external_states['Environment'].nutrients
		self.criticalInitiationMass = self.getDnaCriticalMass(self.nutrientToDoublingTime[current_nutrients])

		##########################################
		# Perform replication initiation process #
		##########################################

		# Calculate how many rounds of replication are occuring on what number of chromosomes
		# information that is used later in the process
		activeDnaPoly = self.activeDnaPoly.molecules()
		activePolymerasePresent = len(activeDnaPoly) > 0

		oriCs = self.oriCs.molecules()
		nChromosomes = self.full_chromosome.total()[0]

		if len(oriCs) == 0 and nChromosomes == 0:
			return

		if activePolymerasePresent:
			replicationRound = activeDnaPoly.attr('replicationRound')

		# Get cell mass
		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg)

		# Calculate mass per origin of replication
		# This is a rearrangement of the equation:
		# 	Cell Mass / Number origin > Critical mass
		# If this is true, initate a round of chromosome replication
		# on every origin of replication
		initiate = False
		massFactor = cellMass / self.criticalInitiationMass
		massPerOrigin = massFactor / len(oriCs)
		if massPerOrigin >= 1.0 and self.chromosomeHalves.total().sum() == 0:
			initiate = True

		# Write data about initiation to listener
		self.writeToListener("ReplicationData", "criticalMassPerOriC", massPerOrigin)
		self.writeToListener("ReplicationData", "criticalInitiationMass", self.criticalInitiationMass.asNumber(units.fg))

		if initiate:
			# If this is triggered initate a round of replication on every origin of replication

			# Calculate number of origins the cell has
			if activePolymerasePresent:
				numOric = 2 * np.unique(replicationRound).size * nChromosomes
			else:
				numOric = 1 * nChromosomes
				replicationRound = np.array([0])

			# Calculate number of new "polymerases" required per origin
			# This is modeled as one "polymerase" on each of the lagging and leading strands.
			# even if there is only a forward and a reverse strand. This was done to make this
			# polymerization process analogous to the transcription and translation elongation processes
			numberOfNewPolymerase = 4 * numOric

			# Initialize 4 "polymerases" per origin
			activeDnaPoly = self.activeDnaPoly.moleculesNew(
				"dnaPolymerase",
				numberOfNewPolymerase
				)

			# Generate a new origin for each old origin
			oriCs = self.oriCs.moleculesNew(
				"originOfReplication",
				numOric
				)

			# Calculate and set attributes of newly created "polymerases"
			sequenceIdx = np.tile(np.array([0,1,2,3], dtype=np.int8), numOric)
			sequenceLength = np.zeros(numberOfNewPolymerase, dtype = np.int8)
			replicationRound = np.ones(numberOfNewPolymerase, dtype=np.int8) * (replicationRound.max() + 1)
			chromosomeIndex = np.zeros(numberOfNewPolymerase, dtype=np.int8)
			chromosomeIndex[numberOfNewPolymerase // nChromosomes:] = 1.

			activeDnaPoly.attrIs(
				sequenceIdx = sequenceIdx,
				sequenceLength = sequenceLength,
				replicationRound = replicationRound,
				chromosomeIndex = chromosomeIndex,
				)


		##########################################
		# Perform replication elongation process #
		##########################################

		# If no "polymerases" are present return
		if len(activeDnaPoly) == 0:
			return

		# Build sequences to polymerize
		dNtpCounts = self.dntps.counts()
		sequenceIdx, sequenceLengths, massDiffDna = activeDnaPoly.attrs(
			'sequenceIdx', 'sequenceLength', 'massDiff_DNA'
			)

		dNtpsUsed = np.zeros_like(dNtpCounts)

		sequences = buildSequences(
			self.sequences,
			sequenceIdx,
			sequenceLengths,
			self._dnaPolymeraseElongationRate()
			)

		# Use polymerize algorithm to quickly calculate the number of elongations
		# each "polymerase" catalyzes
		reactionLimit = dNtpCounts.sum()

		result = polymerize(
			sequences,
			dNtpCounts,
			reactionLimit,
			self.randomState
			)

		sequenceElongations = result.sequenceElongation
		dNtpsUsed = result.monomerUsages
		nElongations = result.nReactions

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

		###########################################
		# Perform replication termination process #
		###########################################

		# Determine if any "polymerases" reached the end of their
		# sequence. If so terminate them and update the attributes of the
		# remaining "polymerases" to reflect the new chromosome structure
		terminalLengths = self.sequenceLengths[sequenceIdx]

		didTerminate = (updatedLengths == terminalLengths)

		terminatedChromosomes = np.bincount(
			sequenceIdx[didTerminate],
			minlength = self.sequences.shape[0]
			)

		if didTerminate.sum():
			# If any forks terminated, update the chromosome index of the next fork up the chromosome
			sequenceIdx, sequenceLength, chromosomeIndex, replicationRound = activeDnaPoly.attrs(
				'sequenceIdx', 'sequenceLength', 'chromosomeIndex', 'replicationRound'
				)

			# Serach for next fork up the chromosome
			# Criteria: Must match the same sequence index, must be current terminating forks's replication
			# round +1, must match terminating fork's chromosome index (so you don't flip the same upstream
			# fork twice).
			sequenceIdxMatch = np.zeros(didTerminate.shape)
			replicationRoundMatch = np.zeros(didTerminate.shape)
			chromosomeIndexMatch = np.zeros(didTerminate.shape)
			for idx in np.where(didTerminate)[0]:
				sequenceIdxMatch = np.logical_or(sequenceIdxMatch, sequenceIdx[idx] == sequenceIdx)
				replicationRoundMatch = np.logical_or(replicationRoundMatch, replicationRound[idx] + 1 == replicationRound)
				chromosomeIndexMatch = np.logical_or(chromosomeIndexMatch, chromosomeIndex[idx] == chromosomeIndex)

			# Potential matches for chromosome index to flip
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

		# Update counts of newly created chromosome halves. These will be "stitched" together
		# in the ChromosomeFormation process
		self.chromosomeHalves.countsInc(terminatedChromosomes)


	def _dnaPolymeraseElongationRate(self):
		# Calculates elongation rate scaled by the time step
		return self.dnaPolyElngRate * self.timeStepSec()

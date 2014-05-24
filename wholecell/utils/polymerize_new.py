"""
polymerize_new.py

Polymerizes sequences based on monomer and energy limitations.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/23/14
"""

from __future__ import division

import numpy as np

PAD_VALUE = -1

# TODO: profile
# TODO: randstream objects

def polymerize(sequences, monomerLimits, reactionLimit):
	# Sanitize inputs
	monomerLimits = monomerLimits.copy().astype(np.int64)
	reactionLimit = np.int64(reactionLimit)

	# Data size information
	nSequences, maxLength = sequences.shape
	nMonomers = monomerLimits.size

	# Static data
	sequenceMonomers = np.rollaxis(np.array([
		sequences == monomerIndex for monomerIndex in xrange(nMonomers)
		]), 0, 3)
	sequenceReactions = (sequences != PAD_VALUE)

	maxElongation = sequenceReactions.sum(1)

	totalBases = sequenceMonomers.sum(0)
	totalReactions = sequenceReactions.sum(0)

	# Running/accumulated values
	activeSequences = sequenceReactions[:, 0]
	cumulativeSequenceMonomers = sequenceMonomers.cumsum(1)
	cumulativeSequenceReactions = sequenceReactions.cumsum(1)

	cumulativeTotalMonomers = totalBases.cumsum(0)
	cumulativeTotalReactions = totalReactions.cumsum(0)

	# Output
	sequenceElongation = np.zeros(nSequences, np.int64)
	monomerUsages = np.zeros(nMonomers, np.int64)
	nReactions = 0

	# Elongate sequences as much as possible
	while activeSequences.any():
		# Find the furthest extent that sequences can elongate without 
		# encountering a limitation

		monomerLimited = [
			np.where(cumulativeTotalMonomers[:, monomerIndex] > monomerLimit)[0]
			for monomerIndex, monomerLimit in enumerate(monomerLimits)
			]

		monomerLimitedAt = np.array([
			limit[0] if limit.size else maxLength
			for limit in monomerLimited
			])

		reactionLimitedAt = maxLength

		reactionLimited = np.where(cumulativeTotalReactions > reactionLimit)[0]

		if reactionLimited.size:
			reactionLimitedAt = reactionLimited[0]

		elongationLimit = np.min([monomerLimitedAt.min(), reactionLimitedAt])

		monomerIsLimiting = (monomerLimitedAt == elongationLimit)
		reactionIsLimited = (reactionLimitedAt == elongationLimit)

		# Elongate up to the limitation

		sequenceElongation[activeSequences] = elongationLimit

		# Update running values

		monomersUsed = cumulativeTotalMonomers[elongationLimit-1]

		monomerUsages += monomersUsed
		monomerLimits -= monomersUsed

		cumulativeSequenceMonomers -= cumulativeSequenceMonomers[:, (elongationLimit-1,), :]
		cumulativeTotalMonomers -= monomersUsed

		reactions = cumulativeTotalReactions[elongationLimit-1]

		nReactions += reactions
		reactionLimit -= reactions

		cumulativeSequenceReactions -= cumulativeSequenceReactions[:, (elongationLimit-1,)]
		cumulativeTotalReactions -= reactions

		# Quit if fully elongated
		if elongationLimit == maxLength:
			break

		# Quit if no remaining monomers
		if monomerLimits.sum() == 0:
			break

		# Quit if no more reaction capacity
		if reactionLimit == 0:
			break

		# Randomly cull sequences that are limiting in monomers
		culledSequences = np.zeros(nSequences, np.bool)

		for monomerIndex, monomerLimit in enumerate(monomerLimits):
			if not monomerIsLimiting[monomerIndex]:
				continue

			activeSequencesWithMonomer = sequenceMonomers[:, elongationLimit, monomerIndex] & activeSequences

			nToCull = activeSequencesWithMonomer.sum() - monomerLimit

			# TODO: pass randstream object
			culledSequenceIndexes = np.random.choice(
				np.where(activeSequencesWithMonomer)[0],
				nToCull,
				replace = False
				)
			
			culledSequences[culledSequenceIndexes] = True

		activeSequences[culledSequences] = False

		# Randomly cull sequences that are limiting in reactions
		if reactionIsLimited:
			activeSequenceWithReaction = sequenceReactions[:, elongationLimit] & activeSequences

			nToCull = activeSequenceWithReaction.sum() - reactionLimit

			# TODO: pass randstream object
			culledSequenceIndexes = np.random.choice(
				np.where(activeSequenceWithReaction)[0],
				nToCull,
				replace = False
				)
			
			culledSequences[culledSequenceIndexes] = True

		activeSequences[culledSequences] = False

		# Update the cumulative totals to include only active sequences

		cumulativeTotalMonomers -= cumulativeSequenceMonomers[culledSequences, :, :].sum(0)
		cumulativeTotalReactions -= cumulativeSequenceReactions[culledSequences, :].sum(0)


	# Restrict elongation by end of sequences

	sequenceElongation = np.fmin(sequenceElongation, maxElongation)

	return sequenceElongation, monomerUsages, nReactions



if __name__ == "__main__":
	import time

	# Contrive a scenario which is similar to real conditions

	nMonomers = 36 # number of distinct aa-tRNAs
	nSequences = 10000 # approximate number of ribosomes
	length = 16 # translation rate
	nTerminating = np.int64(length/300 * nSequences) # estimate for number of ribosomes terminating
	monomerSufficiency = 0.85
	energySufficiency = 0.85

	sequences = np.random.randint(nMonomers, size = (nSequences, length))

	sequenceLengths = length * np.ones(nSequences, np.int64)
	sequenceLengths[np.random.choice(nSequences, nTerminating, replace = False)] = np.random.randint(length, size = nTerminating)

	sequences[np.arange(length) > sequenceLengths[:, np.newaxis]] = PAD_VALUE

	maxReactions = sequenceLengths.sum()

	monomerLimits = maxReactions/nMonomers*np.ones(nMonomers, np.int64)
	reactionLimit = maxReactions

	t = time.time()
	sequenceElongation, monomerUsages, nReactions = polymerize(sequences, monomerLimits, reactionLimit)
	evalTime = time.time() - t

	assert (sequenceElongation <= sequenceLengths+1).all()
	assert (monomerUsages <= monomerLimits).all()
	assert nReactions <= reactionLimit

	print """
Polymerize function report:

For {} sequences of {} different monomers elongating by at most {}:

{:0.1f} ms to evaluate
{} polymerization reactions
{:0.1f} average elongations per sequence
{:0.1%} monomer utilization
{:0.1%} energy utilization
{:0.1%} fully elongated
{:0.1%} completion
""".format(
		nSequences,
		nMonomers,
		length,
		evalTime * 1000,
		nReactions,
		sequenceElongation.mean(),
		monomerUsages.sum()/monomerLimits.sum(),
		nReactions/reactionLimit,
		(sequenceElongation == sequenceLengths).sum()/nSequences,
		sequenceElongation.sum()/sequenceLengths.sum()
		)

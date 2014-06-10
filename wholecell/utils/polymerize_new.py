"""
polymerize_new.py

Polymerizes sequences based on monomer and energy limitations.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/23/14
"""

from __future__ import division

import numpy as np

from _build_sequences import buildSequences, computeMassIncrease

PAD_VALUE = -1

# TODO: consider rewriting as a class to cache sequenceMonomers, which is the
# most expensive operation in the function (alternatively, pass the 
# sequenceMonomers array as an argument instead of sequences)

# TODO: cythonize, since it is starting to look like ~2/3 time is spent in for-loops

def polymerize(sequences, monomerLimits, reactionLimit, randomState):
	# Sanitize inputs
	monomerLimits = monomerLimits.copy().astype(np.int64)
	reactionLimit = np.int64(reactionLimit)

	# Data size information
	nSequences, sequenceLength = sequences.shape
	nMonomers = monomerLimits.size

	# Static data
	sequenceMonomers = np.empty((nMonomers, nSequences, sequenceLength), np.bool)

	for monomerIndex in xrange(nMonomers):
		sequenceMonomers[monomerIndex, ...] = (sequences == monomerIndex)

	sequenceReactions = (sequences != PAD_VALUE)
	sequenceLengths = sequenceReactions.sum(1)

	# Running values
	activeSequencesIndexes = np.arange(nSequences)
	currentStep = 0
	activeSequencesIndexes = activeSequencesIndexes[
		sequenceReactions[:, currentStep]
		]

	totalMonomers = sequenceMonomers.sum(1).cumsum(1)
	totalReactions = sequenceReactions.sum(0).cumsum(0)

	maxElongation = sequenceLength

	# Output
	sequenceElongation = np.zeros(nSequences, np.int64)
	monomerUsages = np.zeros(nMonomers, np.int64)
	nReactions = 0

	# Elongate sequences as much as possible
	while activeSequencesIndexes.size:
		# Find the furthest we can elongate without reaching some limitation
		monomerLimitingExtents = [
			np.where(totalMonomers[monomerIndex, :] > monomerLimit)[0]
			for monomerIndex, monomerLimit in enumerate(monomerLimits)
			]

		monomerLimitedAt = np.array([
			extent[0] if extent.size else maxElongation
			for extent in monomerLimitingExtents
			])

		reactionLimitedAt = maxElongation

		reactionLimitingExtents = np.where(
			totalReactions > reactionLimit
			)[0]

		if reactionLimitingExtents.size:
			reactionLimitedAt = reactionLimitingExtents[0]

		limitingExtent = min(monomerLimitedAt.min(), reactionLimitedAt)

		monomerIsLimiting = (monomerLimitedAt == limitingExtent)
		reactionIsLimiting = (reactionLimitedAt == limitingExtent)

		currentStep += limitingExtent

		# Use resources

		if limitingExtent > 0:
			deltaMonomers = totalMonomers[:, limitingExtent-1]
			deltaReactions = totalReactions[limitingExtent-1]

			monomerLimits -= deltaMonomers
			reactionLimit -= deltaReactions

			monomerUsages += deltaMonomers
			nReactions += deltaReactions

		# Update lengths

		sequenceElongation[activeSequencesIndexes] += limitingExtent

		# Quit if fully elongated

		if limitingExtent == maxElongation:
			break

		# Quit if out of resources

		if not monomerLimits.any():
			break

		if reactionLimit == 0:
			break

		# Cull fully elongated sequences

		sequencesToCull = ~sequenceReactions[activeSequencesIndexes, currentStep]

		# Cull monomer-limiting sequences

		for monomerIndex, monomerLimit in enumerate(monomerLimits):
			if ~monomerIsLimiting[monomerIndex]:
				continue

			sequencesWithMonomer = np.where(
				sequenceMonomers[monomerIndex, activeSequencesIndexes, currentStep]
				)[0]

			nToCull = sequencesWithMonomer.size - monomerLimit

			assert nToCull > 0

			culledIndexes = randomState.choice(
				sequencesWithMonomer,
				nToCull,
				replace = False
				)

			sequencesToCull[culledIndexes] = True

		# Cull reaction-limiting sequences

		if reactionIsLimiting:
			sequencesWithReaction = np.where(
				~sequencesToCull
				)[0]

			nToCull = sequencesWithReaction.size - reactionLimit

			if nToCull > 0:
				culledIndexes = randomState.choice(
					sequencesWithReaction,
					nToCull,
					replace = False
					)

				sequencesToCull[culledIndexes] = True

		# Update running values

		activeSequencesIndexes = activeSequencesIndexes[~sequencesToCull]

		# Quit if there are no more sequences

		if not activeSequencesIndexes.size:
			break

		totalMonomers = sequenceMonomers[:, activeSequencesIndexes, currentStep:].sum(1).cumsum(1)
		totalReactions = sequenceReactions[activeSequencesIndexes, currentStep:].sum(0).cumsum(0)

		maxElongation = sequenceLength - currentStep

	# Clamp sequence lengths up to their max length

	sequenceElongation = np.fmin(sequenceElongation, sequenceLengths)
		
	return sequenceElongation, monomerUsages, nReactions


def _setupExample():
	# Contrive a scenario which is similar to real conditions

	randomState = np.random.RandomState()

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

	monomerLimits = monomerSufficiency * maxReactions/nMonomers*np.ones(nMonomers, np.int64)
	reactionLimit = energySufficiency * maxReactions

	return sequences, monomerLimits, reactionLimit, randomState


def _simpleProfile():
	import time

	np.random.seed(0)

	sequences, monomerLimits, reactionLimit, randomState = _setupExample()

	nSequences, length = sequences.shape
	nMonomers = monomerLimits.size
	sequenceLengths = (sequences != PAD_VALUE).sum(1)

	t = time.time()
	sequenceElongation, monomerUsages, nReactions = polymerize(sequences, monomerLimits, reactionLimit, randomState)
	evalTime = time.time() - t

	assert (sequenceElongation <= sequenceLengths+1).all()
	assert (monomerUsages <= monomerLimits).all()
	assert nReactions <= reactionLimit
	assert nReactions == monomerUsages.sum()

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


def _fullProfile():
	np.random.seed(0)

	sequences, monomerLimits, reactionLimit, randomState = _setupExample()

	# Recipe from https://docs.python.org/2/library/profile.html#module-cProfile

	import cProfile, pstats, StringIO
	pr = cProfile.Profile()
	pr.enable()

	sequenceElongation, monomerUsages, nReactions = polymerize(sequences, monomerLimits, reactionLimit, randomState)
	
	pr.disable()
	s = StringIO.StringIO()
	sortby = 'cumulative'
	ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
	ps.print_stats()
	print s.getvalue()


if __name__ == "__main__":
	_fullProfile()

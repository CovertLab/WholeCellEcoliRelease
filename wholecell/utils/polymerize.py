"""
polymerize.py

Polymerizes sequences based on monomer and energy limitations.

TODO:
- document algorithm/corner cases (should already exist somewhere...)

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/23/14
"""

from __future__ import division

import numpy as np

from _build_sequences import buildSequences, computeMassIncrease
from _fastsums import sum_monomers, sum_monomers_reference_implementation

PAD_VALUE = -1

# Define a no-op @profile decorator if it wasn't loaded by kernprof.
try:
	profile(lambda x: x)
except NameError:
	def profile(function):
		return function

# Run 'kernprof -lv wholecell/tests/utils/profile_polymerize.py' to get a line
# profile of functions decorated @profile.
@profile
def polymerize(sequences, monomerLimits, reactionLimit, randomState):
	# Sanitize inputs
	# import ipdb; ipdb.set_trace()
	monomerLimits = monomerLimits.copy().astype(np.int64)
	reactionLimit = np.int64(reactionLimit)

	# Data size information
	nSequences, sequenceLength = sequences.shape
	nMonomers = monomerLimits.size

	# Static data
	sequenceMonomers = np.empty((nMonomers, nSequences, sequenceLength), dtype = np.bool)

	for monomerIndex in xrange(nMonomers):
		sequenceMonomers[monomerIndex, ...] = (sequences == monomerIndex)

	sequenceReactions = (sequences != PAD_VALUE)
	sequenceLengths = sequenceReactions.sum(axis = 1)

	# Running values
	activeSequencesIndexes = np.arange(nSequences)
	currentStep = 0
	activeSequencesIndexes = activeSequencesIndexes[
		sequenceReactions[:, currentStep]
		]

	#totalMonomers = sequenceMonomers.sum(axis = 1).cumsum(axis = 1)
	totalMonomers = sum_monomers(sequenceMonomers, activeSequencesIndexes, 0)
	totalReactions = sequenceReactions.sum(axis = 0).cumsum(axis = 0)

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

		#totalMonomers = sequenceMonomers[:, activeSequencesIndexes, currentStep:].sum(axis = 1).cumsum(axis = 1)
		#totalMonomers = sum_monomers_reference_implementation(sequenceMonomers, activeSequencesIndexes, currentStep)
		totalMonomers = sum_monomers(sequenceMonomers, activeSequencesIndexes, currentStep)

		totalReactions = sequenceReactions[activeSequencesIndexes, currentStep:].sum(axis = 0).cumsum(axis = 0)

		maxElongation = sequenceLength - currentStep

	# Clamp sequence lengths up to their max length

	sequenceElongation = np.fmin(sequenceElongation, sequenceLengths)
		
	return sequenceElongation, monomerUsages, nReactions

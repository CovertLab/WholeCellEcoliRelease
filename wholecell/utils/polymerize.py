"""
polymerize.py

Polymerizes sequences based on monomer and energy limitations.

Run `kernprof -lv wholecell/tests/utils/profile_polymerize.py` to get a line
profile. It @profile-decorates polymerize().

TODO:
- document algorithm/corner cases (should already exist somewhere...)

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/23/14
"""

from __future__ import absolute_import
from __future__ import division

import numpy as np

from ._build_sequences import buildSequences, computeMassIncrease
from ._fastsums import sum_monomers, sum_monomers_reference_implementation

# Reexport _build_sequences functions. (Declaring this avoids
# "unused import statement" warnings.)
__all__ = ['PAD_VALUE', 'polymerize', 'buildSequences', 'computeMassIncrease']

PAD_VALUE = -1

def polymerize(sequences, monomerLimits, reactionLimit, randomState):
	"""
	Polymerize the given DNA/RNA/protein sequences as far as possible within
	the given limits.

	Parameters:
		sequences: ndarray of integer, shape (num_sequences, num_steps),
			the sequences of needed monomer types, containing PAD_VALUE for all
			steps after sequence completion.
		monomerLimits: ndarray of integer, shape (num_monomers,), the available
			number of each monomer type.
		reactionLimit: max number of reactions (monomers to use); the energy
			limit.
		randomState: random number generator to pick winners in shortages.

	Returns:
		sequenceElongation: ndarray of integer, shape (num_sequences,)
			indicating how far the sequences proceeded,
		monomerUsages: ndarray of integer, shape (num_monomers,) counting how
			many monomers of each type got used,
		nReactions: total number of reactions (monomers used).
	"""
	# Sanitize inputs
	monomerLimits = monomerLimits.astype(np.int64, copy=True)
	reactionLimit = np.int64(reactionLimit)

	# Data size information
	nSequences, sequenceLength = sequences.shape
	nMonomers = monomerLimits.size

	# Static data
	# sequenceMonomers: ndarray of bool, shape
	#     (num_monomers, num_sequences, num_steps), a bitmask of monomer usage.
	sequenceMonomers = np.empty((nMonomers, nSequences, sequenceLength), dtype = np.bool)
	for monomerIndex in xrange(nMonomers):
		sequenceMonomers[monomerIndex, ...] = (sequences == monomerIndex)

	# sequenceReactions: ndarray of bool, shape (num_sequences, num_steps), of
	#     sequence continuation.
	# sequenceLength: ndarray of integer, shape (num_sequences).
	sequenceReactions = (sequences != PAD_VALUE)
	sequenceLengths = sequenceReactions.sum(axis = 1)

	# Running values
	# activeSequencesIndexes: 1D ndarray of integer, the indexes of the
	#     currently active sequences.
	activeSequencesIndexes = np.arange(nSequences)
	currentStep = 0
	activeSequencesIndexes = activeSequencesIndexes[
		sequenceReactions[:, currentStep]
		]

	# totalMonomers: ndarray of integer, shape
	#     (num_monomers, num_steps - currentStep), the count of monomers wanted
	#     in currentStep and beyond.
	# totalReactions: ndarray of integer, shape (num_steps - currentStep,),
	#     the cumulative number of reactions in currentStep to the end for
	#     currently active sequences, ignoring limiting factors.
	#
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
		# Find the furthest we can elongate without reaching some limitation.
		# monomerLimitingExtents: list of ndarray of integer, effective shape
		# (num_monomers, integer), all the currentStep-relative step indexes in
		# currently active sequences beyond the monomer's limit.
		monomerLimitingExtents = [
			np.where(totalMonomers[monomerIndex, :] > monomerLimit)[0]
			for monomerIndex, monomerLimit in enumerate(monomerLimits)
			]

		# monomerLimitedAt: ndarray of integer, shape (num_monomers,), the
		# currentStep-relative step# where each monomer exceeds its limit, else
		# maxElongation.
		monomerLimitedAt = np.array([
			extent[0] if extent.size else maxElongation
			for extent in monomerLimitingExtents
			])

		reactionLimitedAt = maxElongation

		# reactionLimitingExtents: ndarray of integer, shape (integer,), all
		# the currentStep-relative step numbers beyond the reactionLimit.
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
		# sequencesToCull: ndarray of bool, shape (num_active_sequences,),
		# selecting active sequences to cull, initially the ones that finished
		# by currentStep.
		sequencesToCull = ~sequenceReactions[activeSequencesIndexes, currentStep]

		# Cull monomer-limiting sequences
		for monomerIndex, monomerLimit in enumerate(monomerLimits):
			if ~monomerIsLimiting[monomerIndex]:
				continue

			# sequencesWithMonomer: ndarray of integer, shape (integer,), the
			# active sequence indexes that use this monomer in currentStep.
			sequencesWithMonomer = np.where(
				sequenceMonomers[monomerIndex, activeSequencesIndexes, currentStep]
				)[0]

			nToCull = sequencesWithMonomer.size - monomerLimit

			assert nToCull > 0

			# culledIndexes: ndarray of integer, shape (integer,), randomly
			# chosen sequences to cull due to the monomer's limit.
			culledIndexes = randomState.choice(
				sequencesWithMonomer,
				nToCull,
				replace = False
				)

			sequencesToCull[culledIndexes] = True

		# Cull reaction-limiting sequences
		if reactionIsLimiting:
			# sequencesWithReaction: ndarray of integer, shape (integer,), the
			# active sequence indexes not yet ruled out.
			sequencesWithReaction = np.where(
				~sequencesToCull
				)[0]

			nToCull = sequencesWithReaction.size - reactionLimit

			if nToCull > 0:
				# culledIndexes: ndarray of integer, shape (integer,), randomly
				# chosen sequences to cull to uphold reactionLimit.
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

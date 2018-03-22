"""
polymerize.py

Polymerizes sequences based on monomer and energy limitations.

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

# TODO: restore profiling/decoration

class polymerize(object):
	def __init__(self, sequences, monomerLimits, reactionLimit, randomState):
		"""
		Polymerize the given DNA/RNA/protein sequences as far as possible within
		the given limits.

		Args:
			sequences: index[sequence#, step#] sequences of needed monomer types,
				containing PAD_VALUE for all steps after sequence completion
			monomerLimits: count[monomer#] of available monomers
			reactionLimit: scalar max number of reactions (monomers to use)
			randomState: random number generator to pick winners in shortages

		Returns:
			sequenceElongation: length[sequence#] how far the sequences proceeded,
			monomerUsages: count[monomer#] how many monomers got used,
			nReactions: scalar total number of reactions (monomers used)
		"""
		# Sanitize inputs
		monomerLimits = monomerLimits.astype(np.int64, copy=True)
		reactionLimit = np.int64(reactionLimit)

		# Data size information
		nSequences, sequenceLength = sequences.shape
		nMonomers = monomerLimits.size

		# Static data
		# sequenceMonomers: bool[monomer#, sequence#, step#] bitmask of monomer usage
		sequenceMonomers = np.empty((nMonomers, nSequences, sequenceLength), dtype = np.bool)
		for monomerIndex in xrange(nMonomers):
			sequenceMonomers[monomerIndex, ...] = (sequences == monomerIndex)

		# sequenceReactions: bool[sequence#, step#] of sequence continuation
		sequenceReactions = (sequences != PAD_VALUE)
		sequenceLengths = sequenceReactions.sum(axis = 1)

		# Running values
		# activeSequencesIndexes: index[] of sequences that are currently active
		activeSequencesIndexes = np.arange(nSequences)
		currentStep = 0
		activeSequencesIndexes = activeSequencesIndexes[
			sequenceReactions[:, currentStep]
			]

		# totalMonomers:: count:int64[monomer#, step#] of monomers wanted in
		# currentStep and beyond.
		# totalReactions: count:int64[step#] cumulative #reactions
		#
		#totalMonomers = sequenceMonomers.sum(axis = 1).cumsum(axis = 1)
		totalMonomers = sum_monomers(sequenceMonomers, activeSequencesIndexes, 0)
		totalReactions = sequenceReactions.sum(axis = 0).cumsum(axis = 0)

		maxElongation = sequenceLength

		# Output
		self.sequenceElongation = np.zeros(nSequences, np.int64)
		self.monomerUsages = np.zeros(nMonomers, np.int64)
		self.nReactions = 0

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

				self.monomerUsages += deltaMonomers
				self.nReactions += deltaReactions

			# Update lengths
			self.sequenceElongation[activeSequencesIndexes] += limitingExtent

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
		self.sequenceElongation = np.fmin(self.sequenceElongation, sequenceLengths)

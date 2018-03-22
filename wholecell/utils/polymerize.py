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
__all__ = ['polymerize', 'buildSequences', 'computeMassIncrease']


# TODO (John): Restore line profiler decorator?

class polymerize(object): # Class name is lowercase because interface is function-like
	PAD_VALUE = -1
	def __init__(self, sequences, monomerLimits, reactionLimit, randomState):
		# TODO (John): should the docstring be under __init__ or the class?
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

		# Gather inputs
		self._sequences = sequences
		self._monomerLimits = monomerLimits
		self._reactionLimit = reactionLimit
		self._randomState = randomState

		# Prepare for iteration
		self._sanitize_inputs()
		self._gather_input_dimensions()
		self._gather_sequence_data()
		# self._prepare_running_values()

		# Running values
		# activeSequencesIndexes: index[] of sequences that are currently active
		activeSequencesIndexes = np.arange(self._nSequences)
		currentStep = 0
		activeSequencesIndexes = activeSequencesIndexes[
			self._sequenceReactions[:, currentStep]
			]

		# totalMonomers:: count:int64[monomer#, step#] of monomers wanted in
		# currentStep and beyond.
		# totalReactions: count:int64[step#] cumulative #reactions
		#
		#totalMonomers = self._sequenceMonomers.sum(axis = 1).cumsum(axis = 1)
		totalMonomers = sum_monomers(self._sequenceMonomers, activeSequencesIndexes, 0)
		totalReactions = self._sequenceReactions.sum(axis = 0).cumsum(axis = 0)

		maxElongation = self._sequenceLength

		# Output
		self.sequenceElongation = np.zeros(self._nSequences, np.int64)
		self.monomerUsages = np.zeros(self._nMonomers, np.int64)
		self.nReactions = 0

		# Elongate sequences as much as possible
		while activeSequencesIndexes.size:
			# Find the furthest we can elongate without reaching some limitation
			monomerLimitingExtents = [
				np.where(totalMonomers[monomerIndex, :] > monomerLimit)[0]
				for monomerIndex, monomerLimit in enumerate(self._monomerLimits)
				]

			monomerLimitedAt = np.array([
				extent[0] if extent.size else maxElongation
				for extent in monomerLimitingExtents
				])

			reactionLimitedAt = maxElongation

			reactionLimitingExtents = np.where(
				totalReactions > self._reactionLimit
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

				self._monomerLimits -= deltaMonomers
				self._reactionLimit -= deltaReactions

				self.monomerUsages += deltaMonomers
				self.nReactions += deltaReactions

			# Update lengths
			self.sequenceElongation[activeSequencesIndexes] += limitingExtent

			# Quit if fully elongated
			if limitingExtent == maxElongation:
				break

			# Quit if out of resources
			if not self._monomerLimits.any():
				break

			if self._reactionLimit == 0:
				break

			# Cull fully elongated sequences
			sequencesToCull = ~self._sequenceReactions[activeSequencesIndexes, currentStep]

			# Cull monomer-limiting sequences
			for monomerIndex, monomerLimit in enumerate(self._monomerLimits):
				if ~monomerIsLimiting[monomerIndex]:
					continue

				sequencesWithMonomer = np.where(
					self._sequenceMonomers[monomerIndex, activeSequencesIndexes, currentStep]
					)[0]

				nToCull = sequencesWithMonomer.size - monomerLimit

				assert nToCull > 0

				culledIndexes = self._randomState.choice(
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

				nToCull = sequencesWithReaction.size - self._reactionLimit

				if nToCull > 0:
					culledIndexes = self._randomState.choice(
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

			#totalMonomers = self._sequenceMonomers[:, activeSequencesIndexes, currentStep:].sum(axis = 1).cumsum(axis = 1)
			#totalMonomers = sum_monomers_reference_implementation(self._sequenceMonomers, activeSequencesIndexes, currentStep)
			totalMonomers = sum_monomers(self._sequenceMonomers, activeSequencesIndexes, currentStep)

			totalReactions = self._sequenceReactions[activeSequencesIndexes, currentStep:].sum(axis = 0).cumsum(axis = 0)

			maxElongation = self._sequenceLength - currentStep

		# Clamp sequence lengths up to their max length
		self.sequenceElongation = np.fmin(self.sequenceElongation, self._sequenceLengths)

	# __init__ subroutines
	# Several of these assign new attributes outside of __init__'s context, but
	# these functions should never be called more than once or outside __init__
	# as they are just part of extended __init__ operations.

	def _sanitize_inputs(self):
		'''
		Enforce array typing, and copy input arrays to prevent side-effects.
		'''
		self._monomerLimits = self._monomerLimits.astype(np.int64, copy = True)
		self._reactionLimit = np.int64(self._reactionLimit)

	def _gather_input_dimensions(self):
		'''
		Collect information about the size of the inputs.
		'''

		(self._nSequences, self._sequenceLength) = self._sequences.shape
		self._nMonomers = self._monomerLimits.size

	def _gather_sequence_data(self):
		'''
		Collect static data about the input sequences.
		'''

		# self._sequenceMonomers: bool[monomer#, sequence#, step#] bitmask of monomer usage
		self._sequenceMonomers = np.empty(
			(self._nMonomers, self._nSequences, self._sequenceLength),
			dtype = np.bool
			)
		for monomerIndex in xrange(self._nMonomers):
			self._sequenceMonomers[monomerIndex, ...] = (self._sequences == monomerIndex)

		# self._sequenceReactions: bool[sequence#, step#] of sequence continuation
		self._sequenceReactions = (self._sequences != self.PAD_VALUE)
		self._sequenceLengths = self._sequenceReactions.sum(axis = 1)

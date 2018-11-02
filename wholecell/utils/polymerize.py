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
__all__ = ['polymerize', 'buildSequences', 'computeMassIncrease']

class polymerize(object): # Class name is lowercase because interface is function-like
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

	PAD_VALUE = -1

	def __init__(self, sequences, monomerLimits, reactionLimit, randomState):
		# Gather inputs
		self._sequences = sequences
		self._monomerLimits = monomerLimits
		self._reactionLimit = reactionLimit
		self._randomState = randomState

		# Prepare for iteration
		self._setup()

		# Elongate sequences as much as possible
		self._elongate()

		# Clean up
		self._finalize()

	# __init__ subroutines
	# Several of these assign new attributes outside of __init__'s immediate
	# context; however, they should only ever be called by __init__.

	# Setup subroutines

	def _setup(self):
		'''
		Extended initialization procedures.
		'''
		self._sanitize_inputs()
		self._gather_input_dimensions()
		self._gather_sequence_data()
		self._prepare_running_values()
		self._prepare_outputs()

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

		# sequenceMonomers: ndarray of bool, shape
		#     (num_monomers, num_sequences, num_steps), a bitmask of monomer usage.
		self._sequenceMonomers = np.empty(
			(self._nMonomers, self._nSequences, self._sequenceLength),
			dtype = np.bool
			)
		for monomerIndex in xrange(self._nMonomers):
			self._sequenceMonomers[monomerIndex, ...] = (
				self._sequences == monomerIndex
				)

		# sequenceReactions: ndarray of bool, shape (num_sequences, num_steps), of
		#     sequence continuation.
		# sequenceLength: ndarray of integer, shape (num_sequences).
		self._sequenceReactions = (self._sequences != self.PAD_VALUE)
		self._sequenceLengths = self._sequenceReactions.sum(axis = 1)

	def _prepare_running_values(self):
		'''
		Sets up the variables that will change throughout iteration, including
		both intermediate calculations and outputs.
		'''

		# activeSequencesIndexes: 1D ndarray of integer, the indexes of the
		#     currently active sequences.
		self._activeSequencesIndexes = np.arange(self._nSequences)
		self._currentStep = 0
		self._activeSequencesIndexes = self._activeSequencesIndexes[
			self._sequenceReactions[:, self._currentStep]
			]

		self._update_elongation_resource_demands()
		self._monomerHistory = np.empty((self._maxElongation, self._nMonomers), np.int64)

		# Empty placeholders - will be filled in during trivial elongation,
		# then inspected during nontrivial (resource-limited) elongation
		self._monomerIsLimiting = np.empty(self._nMonomers, np.bool)
		self._reactionIsLimiting = None

	def _prepare_outputs(self):
		'''
		Running values that ultimately compose the output of the 'polymerize'
		operation.
		'''

		self.sequenceElongation = np.zeros(self._nSequences, np.int64)
		self.monomerUsages = np.zeros(self._nMonomers, np.int64)
		self.nReactions = 0

	# Iteration subroutines

	def _elongate(self):
		'''
		Iteratively elongates sequences up to resource limits.
		'''

		while True:
			# Perform trivial elongations
			fully_elongated = self._elongate_to_limit()

			monomer_limited = (self._monomerLimits == 0).all()
			reaction_limited = (self._reactionLimit == 0)

			# Quit if finished or out of resources
			if fully_elongated or monomer_limited or reaction_limited:
				break

			# Perform nontrivial (resource-limited) elongations, and cull
			# sequences that can no longer be elongated

			self._finalize_resource_limited_elongations()

			# Quit if there are no more sequences
			if not self._activeSequencesIndexes.size:
				break

			# Otherwise, update running values
			self._update_elongation_resource_demands()

	def _elongate_to_limit(self):
		'''
		Elongate as far as possible without hitting any resource limitations.
		'''

		projectionIndex = 0
		notLimited = True
		monomerProjection = np.zeros(self._nMonomers, np.int64)

		# advance one step at a time until a sequence is limited
		while notLimited and projectionIndex < self._maxElongation:
			index = self._currentStep + projectionIndex
			monomerStep = sum_monomers(self._sequenceMonomers[:, :, index], self._activeSequencesIndexes)
			self._monomerHistory[projectionIndex] = monomerProjection + monomerStep
			self._monomerIsLimiting = self._monomerHistory[projectionIndex] > self._monomerLimits
			if self._monomerIsLimiting.any():
				notLimited = False
			else:
				monomerProjection += monomerStep
				projectionIndex += 1

		reactionLimitedAt = self._maxElongation

		# reactionLimitingExtents: ndarray of integer, shape (integer,), all
		# the currentStep-relative step numbers beyond the reactionLimit.
		reactionLimitingExtents = np.where(
			self._totalReactions > self._reactionLimit
			)[0]

		if reactionLimitingExtents.size:
			reactionLimitedAt = reactionLimitingExtents[0]

		limitingExtent = min(projectionIndex, reactionLimitedAt)

		self._reactionIsLimiting = (reactionLimitedAt == limitingExtent)

		self._currentStep += limitingExtent

		# Use resources
		if limitingExtent > 0:
			deltaMonomers = self._monomerHistory[limitingExtent-1]
			deltaReactions = self._totalReactions[limitingExtent-1]

			self._monomerLimits -= deltaMonomers
			self._reactionLimit -= deltaReactions

			self.monomerUsages += deltaMonomers
			self.nReactions += deltaReactions

		# Update lengths
		self.sequenceElongation[self._activeSequencesIndexes] += limitingExtent

		# Determine whether we are finished elongating
		# TODO (John): see if we can determine this outside this context,
		# consequently removing the need to "return" anything
		fully_elongated = (limitingExtent == self._maxElongation)

		return fully_elongated

	def _finalize_resource_limited_elongations(self):
		# Find fully elongated sequences

		# sequencesToCull: ndarray of bool, shape (num_active_sequences,),
		# selecting active sequences to cull, initially the ones that finished
		# by currentStep.
		sequencesToCull = ~self._sequenceReactions[self._activeSequencesIndexes, self._currentStep]

		# Find and finalize monomer-limiting sequences
		for monomerIndex, monomerLimit in enumerate(self._monomerLimits):
			if ~self._monomerIsLimiting[monomerIndex]:
				continue

			# sequencesWithMonomer: ndarray of integer, shape (integer,), the
			# active sequence indexes that use this monomer in currentStep.
			sequencesWithMonomer = np.where(
				self._sequenceMonomers[monomerIndex, self._activeSequencesIndexes, self._currentStep]
				)[0]

			nToCull = sequencesWithMonomer.size - monomerLimit

			assert nToCull > 0

			# culledIndexes: ndarray of integer, shape (integer,), randomly
			# chosen sequences to cull due to the monomer's limit.
			culledIndexes = self._randomState.choice(
				sequencesWithMonomer,
				nToCull,
				replace = False
				)

			sequencesToCull[culledIndexes] = True

		# Find and finalize reaction-limiting sequences
		if self._reactionIsLimiting:
			# sequencesWithReaction: ndarray of integer, shape (integer,), the
			# active sequence indexes not yet ruled out.
			sequencesWithReaction = np.where(
				~sequencesToCull
				)[0]

			nToCull = sequencesWithReaction.size - self._reactionLimit

			if nToCull > 0:
				# culledIndexes: ndarray of integer, shape (integer,), randomly
				# chosen sequences to cull to uphold reactionLimit.
				culledIndexes = self._randomState.choice(
					sequencesWithReaction,
					nToCull,
					replace = False
					)

				sequencesToCull[culledIndexes] = True

		# Cull sequences
		self._activeSequencesIndexes = self._activeSequencesIndexes[~sequencesToCull]

	def _update_elongation_resource_demands(self):
		'''
		After updating the active sequences (initialization and culling),
		recalculate resource demands for the remaining steps given what
		sequences remain.
		'''

		# totalReactions: ndarray of integer, shape (num_steps - currentStep,),
		#     the cumulative number of reactions in currentStep to the end for
		#     currently active sequences, ignoring limiting factors.
		self._totalReactions = self._sequenceReactions[self._activeSequencesIndexes, self._currentStep:].sum(axis = 0).cumsum(axis = 0)

		self._maxElongation = self._sequenceLength - self._currentStep

	# Finalization subroutines

	def _finalize(self):
		'''
		Clean up iteration results.
		'''

		self._clamp_elongation_to_sequence_length()

	def _clamp_elongation_to_sequence_length(self):
		'''
		A post-iteration clean-up operation.  Restricts the elongation of a
		sequence to at most its total (unpadded) length.

		TODO: explain why we do this here instead of during each iteration
		'''

		self.sequenceElongation = np.fmin(
			self.sequenceElongation,
			self._sequenceLengths
			)

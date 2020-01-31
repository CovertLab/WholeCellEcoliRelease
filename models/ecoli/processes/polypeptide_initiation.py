#!/usr/bin/env python
"""
PolypeptideInitiation

Polypeptide initiation sub-model.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/30/14
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils import units
from wholecell.utils.fitting import normalize


class PolypeptideInitiation(wholecell.processes.process.Process):
	""" PolypeptideInitiation """

	_name = "PolypeptideInitiation"

	def __init__(self):
		super(PolypeptideInitiation, self).__init__()

	def initialize(self, sim, sim_data):
		super(PolypeptideInitiation, self).initialize(sim, sim_data)

		# Load parameters
		self.proteinLengths = sim_data.process.translation.monomerData["length"].asNumber()
		self.translationEfficiencies = normalize(
			sim_data.process.translation.translationEfficienciesByMonomer)
		self.fracActiveRibosomeDict = sim_data.process.translation.ribosomeFractionActiveDict
		self.ribosomeElongationRateDict = sim_data.process.translation.ribosomeElongationRateDict
		self.variable_elongation = sim._variable_elongation_translation
		self.make_elongation_rates = sim_data.process.translation.make_elongation_rates

		# Get indexes from proteins to transcription units
		self.protein_index_to_TU_index = sim_data.relation.rnaIndexToMonomerMapping

		# Build matrix to convert transcription unit counts to mRNA counts
		all_TU_ids = sim_data.process.transcription.rnaData['id']
		all_mRNA_ids = sim_data.process.translation.monomerData['rnaId']
		self.n_TUs = len(all_TU_ids)
		self.n_mRNAs = len(all_mRNA_ids)

		self.TU_counts_to_mRNA_counts = np.zeros(
			(self.n_mRNAs, self.n_TUs), dtype=np.float64)

		TU_id_to_index = {TU_id: i for i, TU_id in enumerate(all_TU_ids)}
		for i, mRNA_id in enumerate(all_mRNA_ids):
			self.TU_counts_to_mRNA_counts[i, TU_id_to_index[mRNA_id]] = 1

		# Determine changes from parameter shuffling variant
		if (hasattr(sim_data.process.translation, "translationEfficienciesShuffleIdxs")
				and sim_data.process.translation.translationEfficienciesShuffleIdxs is not None):
			shuffleIdxs = sim_data.process.translation.translationEfficienciesShuffleIdxs
			self.translationEfficiencies = self.translationEfficiencies[shuffleIdxs]

		# Create view on to active 70S ribosomes
		self.active_ribosomes = self.uniqueMoleculesView('active_ribosome')

		# Create views onto bulk 30S and 50S ribosomal subunits
		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeIds.s30_fullComplex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeIds.s50_fullComplex)

		# Create view onto RNAs
		self.RNAs = self.uniqueMoleculesView('RNA')

	def calculateRequest(self):
		current_media_id = self._external_states['Environment'].current_media_id

		self.ribosome30S.requestAll()
		self.ribosome50S.requestAll()

		self.fracActiveRibosome = self.fracActiveRibosomeDict[current_media_id]

		# Read ribosome elongation rate from last timestep
		self.ribosomeElongationRate = self.readFromListener(
			"RibosomeData", "effectiveElongationRate")
		# If the ribosome elongation rate is zero (which is always the case for
		# the first timestep), set ribosome elongation rate to the one in
		# dictionary
		if self.ribosomeElongationRate == 0:
			self.ribosomeElongationRate = self.ribosomeElongationRateDict[
				current_media_id].asNumber(units.aa / units.s)
		self.elongation_rates = self.make_elongation_rates(
			self.randomState,
			self.ribosomeElongationRate,
			1,  # want elongation rate, not lengths adjusted for time step
			self.variable_elongation)

		# Ensure rates are never zero
		self.elongation_rates = np.fmax(self.elongation_rates, 1)

	def evolveState(self):
		# Calculate number of ribosomes that could potentially be initialized
		# based on counts of free 30S and 50S subunits
		inactiveRibosomeCount = np.min([
			self.ribosome30S.count().sum(),
			self.ribosome50S.count().sum(),
			])

		# Get attributes of active (translatable) mRNAs
		TU_index_RNAs, can_translate, unique_index_RNAs = self.RNAs.attrs(
			'TU_index', 'can_translate', 'unique_index')
		TU_index_active_mRNAs = TU_index_RNAs[can_translate]
		unique_index_active_mRNAs = unique_index_RNAs[can_translate]

		# Get counts of each type of active mRNA
		TU_counts = np.bincount(TU_index_active_mRNAs, minlength=self.n_TUs)
		mRNA_counts = self.TU_counts_to_mRNA_counts.dot(TU_counts)

		# Calculate initiation probabilities for ribosomes based on mRNA counts
		# and associated mRNA translational efficiencies
		proteinInitProb = normalize(
			mRNA_counts * self.translationEfficiencies
			)

		# Calculate actual number of ribosomes that should be activated based
		# on probabilities
		self.activationProb = self._calculateActivationProb(
			self.fracActiveRibosome,
			self.proteinLengths,
			self.elongation_rates,
			proteinInitProb,
			self.timeStepSec())

		n_ribosomes_to_activate = np.int64(self.activationProb * inactiveRibosomeCount)

		if n_ribosomes_to_activate == 0:
			return

		# Sample multinomial distribution to determine which mRNAs have full
		# 70S ribosomes initialized on them
		n_new_proteins = self.randomState.multinomial(
			n_ribosomes_to_activate,
			proteinInitProb
			)

		# Build attributes for active ribosomes.
		# Each ribosome is assigned a protein index for the protein that
		# corresponds to the polypeptide it will polymerize. This is done in
		# blocks of protein ids for efficiency.
		protein_indexes = np.empty(n_ribosomes_to_activate, np.int64)
		mRNA_indexes = np.empty(n_ribosomes_to_activate, np.int64)
		nonzeroCount = (n_new_proteins > 0)
		start_index = 0

		for protein_index, counts in zip(
				np.arange(n_new_proteins.size)[nonzeroCount],
				n_new_proteins[nonzeroCount]):
			# Set protein index
			protein_indexes[start_index:start_index+counts] = protein_index

			# Get mask for active mRNA molecules that produce this protein
			mask = (TU_index_active_mRNAs == self.protein_index_to_TU_index[protein_index])
			n_mRNAs = np.count_nonzero(mask)

			# Distribute ribosomes among these mRNAs
			n_ribosomes_per_RNA = self.randomState.multinomial(
				counts, np.full(n_mRNAs, 1./n_mRNAs))

			# Get unique indexes of each mRNA
			mRNA_indexes[start_index:start_index + counts] = np.repeat(
				unique_index_active_mRNAs[mask], n_ribosomes_per_RNA)

			start_index += counts

		# Create active 70S ribosomes and assign their attributes
		self.active_ribosomes.moleculesNew(
			n_ribosomes_to_activate,
			protein_index=protein_indexes,
			peptide_length=np.zeros(n_ribosomes_to_activate, dtype=np.int64),
			mRNA_index=mRNA_indexes,
			pos_on_mRNA=np.zeros(n_ribosomes_to_activate, dtype=np.int64)
			)

		# Decrement free 30S and 70S ribosomal subunit counts
		self.ribosome30S.countDec(n_new_proteins.sum())
		self.ribosome50S.countDec(n_new_proteins.sum())

		# Write number of initialized ribosomes to listener
		self.writeToListener("RibosomeData", "didInitialize", n_new_proteins.sum())
		self.writeToListener("RibosomeData", "probTranslationPerTranscript", proteinInitProb)

	def _calculateActivationProb(self, fracActiveRibosome, proteinLengths, ribosomeElongationRates, proteinInitProb, timeStepSec):
		"""
		Calculates the expected ribosome termination rate based on the ribosome
		elongation rate
		Params:
			- allTranslationTimes: Vector of times required to translate each
			protein
			- allTranslationTimestepCounts: Vector of numbers of timesteps
			required to translate each protein
			- averageTranslationTimeStepCounts: Average number of timesteps
			required to translate a protein, weighted by initiation
			probabilities
			- expectedTerminationRate: Average number of terminations in one
			timestep for one protein
		"""
		allTranslationTimes = 1. / ribosomeElongationRates * proteinLengths
		allTranslationTimestepCounts = np.ceil(allTranslationTimes / timeStepSec)
		averageTranslationTimestepCounts = np.dot(allTranslationTimestepCounts, proteinInitProb)
		expectedTerminationRate = 1.0 / averageTranslationTimestepCounts

		# Modify given fraction of active ribosomes to take into account early
		# terminations in between timesteps
		# allFractionTimeInactive: Vector of probabilities an "active" ribosome
		# 	will in effect be "inactive" because it has terminated during a
		# 	timestep
		# averageFractionTimeInactive: Average probability of an "active"
		# 	ribosome being in effect "inactive", weighted by initiation
		#	probabilities
		# effectiveFracActiveRnap: New higher "goal" for fraction of active
		# 	ribosomes, considering that the "effective" fraction is lower than
		# 	what the listener sees
		allFractionTimeInactive = 1 - allTranslationTimes / timeStepSec / allTranslationTimestepCounts
		averageFractionTimeInactive = np.dot(allFractionTimeInactive, proteinInitProb)
		effectiveFracActiveRibosome = fracActiveRibosome * 1 / (1 - averageFractionTimeInactive)

		# Return activation probability that will balance out the expected
		# termination rate
		activationProb = effectiveFracActiveRibosome * expectedTerminationRate / (1 - effectiveFracActiveRibosome)

		# The upper bound for the activation probability is temporarily set to
		# 1.0 to prevent negative molecule counts. This will lower the fraction
		# of active ribosomes for timesteps longer than roughly 1.8s.
		if activationProb >= 1.0:
			activationProb = 1

		return activationProb

	def isTimeStepShortEnough(self, inputTimeStep, timeStepSafetyFraction):
		# Return false if timestep is longer than 1.7s. Timesteps longer than
		# 1.8s will lower the fraction of activated ribosomes.
		if inputTimeStep > 1.7:
			return False
		else:
			return True

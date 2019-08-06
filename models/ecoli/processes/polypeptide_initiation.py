#!/usr/bin/env python

"""
PolypeptideInitiation

Polypeptide initiation sub-model.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/30/14
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils import units
from wholecell.utils.fitting import normalize

import itertools

class PolypeptideInitiation(wholecell.processes.process.Process):
	""" PolypeptideInitiation """

	_name = "PolypeptideInitiation"

	def __init__(self):
		super(PolypeptideInitiation, self).__init__()

	def initialize(self, sim, sim_data):
		super(PolypeptideInitiation, self).initialize(sim, sim_data)

		# Load parameters
		mrnaIds = sim_data.process.translation.monomerData["rnaId"]
		self.proteinLengths = sim_data.process.translation.monomerData["length"].asNumber()
		self.translationEfficiencies = normalize(sim_data.process.translation.translationEfficienciesByMonomer)
		self.fracActiveRibosomeDict = sim_data.process.translation.ribosomeFractionActiveDict
		self.ribosomeElongationRateDict = sim_data.process.translation.ribosomeElongationRateDict
		self.variable_elongation = sim._variable_elongation_translation
		self.make_elongation_rates = sim_data.process.translation.make_elongation_rates

		# Determine changes from parameter shuffling variant
		shuffleIdxs = None
		if hasattr(sim_data.process.translation, "translationEfficienciesShuffleIdxs") and sim_data.process.translation.translationEfficienciesShuffleIdxs != None:
			shuffleIdxs = sim_data.process.translation.translationEfficienciesShuffleIdxs
			self.translationEfficiencies = self.translationEfficiencies[shuffleIdxs]

		# Create view on to active 70S ribosomes
		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')

		# Create views onto bulk 30S and 50S ribosomal subunits
		self.ribosome30S = self.bulkMoleculeView(sim_data.moleculeIds.s30_fullComplex)
		self.ribosome50S = self.bulkMoleculeView(sim_data.moleculeIds.s50_fullComplex)

		# Create view onto bulk mRNAs
		self.mRnas = self.bulkMoleculesView(mrnaIds)

	def calculateRequest(self):
		current_media_id = self._external_states['Environment'].current_media_id

		self.ribosome30S.requestAll()
		self.ribosome50S.requestAll()
		self.mRnas.requestAll()

		self.fracActiveRibosome = self.fracActiveRibosomeDict[current_media_id]

		# Read ribosome elongation rate from last timestep
		self.ribosomeElongationRate = self.readFromListener("RibosomeData", "effectiveElongationRate")
		# If the ribosome elongation rate is zero (which is always the case for the first timestep), set ribosome elongation rate to one in dictionary
		if self.ribosomeElongationRate == 0:
			self.ribosomeElongationRate = self.ribosomeElongationRateDict[current_media_id].asNumber(units.aa / units.s)
		self.elongation_rates = self.make_elongation_rates(
			self.randomState,
			self.ribosomeElongationRate,
			self.timeStepSec(),
			self.variable_elongation)


	def evolveState(self):
		# Calculate number of ribosomes that could potentially be initalized based on
		# counts of free 30S and 50S subunits
		inactiveRibosomeCount = np.min([
			self.ribosome30S.count().sum(),
			self.ribosome50S.count().sum(),
			])

		# Calculate initiation probabilities for ribosomes based on mRNA counts and associated
		# mRNA translational efficiencies
		proteinInitProb = normalize(
			self.mRnas.counts() * self.translationEfficiencies
			)


		# Calculate actual number of ribosomes that should be activated based on probabilities
		self.activationProb = self._calculateActivationProb(
			self.fracActiveRibosome,
			self.proteinLengths,
			self.elongation_rates,
			proteinInitProb,
			self.timeStepSec())

		ribosomeToActivate = np.int64(self.activationProb * inactiveRibosomeCount)

		if ribosomeToActivate == 0:
			return

		# Sample multinomial distribution to determine which mRNAs have full 70S
		# ribosomes initalized on them
		nNewProteins = self.randomState.multinomial(
			ribosomeToActivate,
			proteinInitProb
			)

		# Each ribosome is assigned a protein index for the protein that corresponds to the
		# polypeptide it will polymerize. This is done in blocks of protein ids for efficiency.
		proteinIndexes = np.empty(ribosomeToActivate, np.int64)
		nonzeroCount = (nNewProteins > 0)
		startIndex = 0
		for proteinIndex, counts in itertools.izip(
				np.arange(nNewProteins.size)[nonzeroCount],
				nNewProteins[nonzeroCount],
				):

			proteinIndexes[startIndex:startIndex+counts] = proteinIndex
			startIndex += counts

		# Create active 70S ribosomes and assign their protein indexes calculated above
		self.activeRibosomes.moleculesNew(
			ribosomeToActivate,
			proteinIndex = proteinIndexes
			)

		# Decrement free 30S and 70S ribosomal subunit counts
		self.ribosome30S.countDec(nNewProteins.sum())
		self.ribosome50S.countDec(nNewProteins.sum())

		# Write number of initalized ribosomes to listener
		self.writeToListener("RibosomeData", "didInitialize", nNewProteins.sum())
		self.writeToListener("RibosomeData", "probTranslationPerTranscript", proteinInitProb)

	def _calculateActivationProb(self, fracActiveRibosome, proteinLengths, ribosomeElongationRates, proteinInitProb, timeStepSec):
		# Calculate expected ribosome termination rate based on ribosome elongation rate
		# allTranslationTimes: Vector of times required to translate each protein
		# allTranslationTimestepCounts: Vector of numbers of timesteps required to translate each protein
		# averageTranslationTimeStepCounts: Average number of timesteps required to translate a protein, weighted by initiation probabilities
		# expectedTerminationRate: Average number of terminations in one timestep for one protein
		allTranslationTimes = 1. / ribosomeElongationRates * proteinLengths
		allTranslationTimestepCounts = np.ceil(allTranslationTimes / (timeStepSec * 1.0))
		averageTranslationTimestepCounts = np.dot(allTranslationTimestepCounts, proteinInitProb)
		expectedTerminationRate = 1.0 / averageTranslationTimestepCounts

		# Modify given fraction of active ribosomes to take into account early terminations in between timesteps
		# allFractionTimeInactive: Vector of probabilities an "active" ribosome will in effect be "inactive" because it has terminated during a timestep
		# averageFractionTimeInactive: Average probability of an "active" ribosome being in effect "inactive", weighted by initiation probabilities
		# effectiveFracActiveRnap: New higher "goal" for fraction of active ribosomes, considering that the "effective" fraction is lower than what the listener sees
		allFractionTimeInactive = 1 - allTranslationTimes / (timeStepSec * 1.0) / allTranslationTimestepCounts
		averageFractionTimeInactive = np.dot(allFractionTimeInactive, proteinInitProb)
		effectiveFracActiveRibosome = fracActiveRibosome * 1 / (1 - averageFractionTimeInactive)

		# Return activation probability that will balance out the expected termination rate
		activationProb = effectiveFracActiveRibosome * expectedTerminationRate / (1 - effectiveFracActiveRibosome)

		# The upper bound for the activation probability is temporarily set to 1.0 to prevent negative molecule counts. This will
		# lower the fraction of active ribosomes for timesteps longer than roughly 1.8s.
		if activationProb >= 1.0:
			activationProb = 1

		return activationProb

	def isTimeStepShortEnough(self, inputTimeStep, timeStepSafetyFraction):
		# Return false if timestep is longer than 1.7s. Timesteps longer than 1.8s will lower the fraction of activated ribosomes
		if inputTimeStep > 1.7:
			return False
		else:
			return True

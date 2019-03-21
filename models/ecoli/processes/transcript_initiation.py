#!/usr/bin/env python

"""
TranscriptInitiation

Transcription initiation sub-model.

TODO:
- use transcription units instead of single genes
- match sigma factors to promoters

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/26/14
"""

from __future__ import division

import numpy as np
import scipy.sparse

import wholecell.processes.process
from wholecell.utils import units

from itertools import izip

class TranscriptInitiation(wholecell.processes.process.Process):
	""" TranscriptInitiation """

	_name = "TranscriptInitiation"


	# Constructor
	def __init__(self):
		super(TranscriptInitiation, self).__init__()

	def initialize(self, sim, sim_data):
		super(TranscriptInitiation, self).initialize(sim, sim_data)

		# Load parameters
		self.fracActiveRnapDict = sim_data.process.transcription.rnapFractionActiveDict
		self.rnaLengths = sim_data.process.transcription.rnaData["length"]
		self.rnaPolymeraseElongationRateDict = sim_data.process.transcription.rnaPolymeraseElongationRateDict

		# Initialize matrices used to calculate synthesis probabilities
		self.basal_prob = sim_data.process.transcription_regulation.basal_prob
		self.n_TUs = len(self.basal_prob)
		delta_prob = sim_data.process.transcription_regulation.delta_prob
		self.delta_prob_matrix = scipy.sparse.csr_matrix(
			(delta_prob['deltaV'],
			(delta_prob['deltaI'], delta_prob['deltaJ'])),
			shape=delta_prob['shape']
			).toarray()

		self.maxRibosomeElongationRate = float(
			sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s))

		# Determine changes from genetic perturbations
		self.genetic_perturbations = {}
		perturbations = getattr(sim_data, "genetic_perturbations", {})

		if len(perturbations) > 0:
			probability_indexes = [
				(index, sim_data.genetic_perturbations[rna_data['id']])
					for index, rna_data in enumerate(sim_data.process.transcription.rnaData)
					if rna_data['id'] in sim_data.genetic_perturbations]

			self.genetic_perturbations = {
				'fixedRnaIdxs': map(lambda pair: pair[0], probability_indexes),
				'fixedSynthProbs': map(lambda pair: pair[1], probability_indexes)
				}

		# If initiationShuffleIdxs does not exist, set value to None
		self.shuffleIdxs = getattr(
			sim_data.process.transcription, "initiationShuffleIdxs", None)

		# Views
		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')
		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")
		self.full_chromosomes = self.uniqueMoleculesView('fullChromosome')
		self.promoters = self.uniqueMoleculesView('promoter')

		# ID Groups
		self.idx_16Srrna = np.where(sim_data.process.transcription.rnaData['isRRna16S'])[0]
		self.idx_23Srrna = np.where(sim_data.process.transcription.rnaData['isRRna23S'])[0]
		self.idx_5Srrna = np.where(sim_data.process.transcription.rnaData['isRRna5S'])[0]
		self.idx_rrna = np.where(sim_data.process.transcription.rnaData['isRRna'])[0]
		self.idx_mrna = np.where(sim_data.process.transcription.rnaData["isMRna"])[0]
		self.idx_trna = np.where(sim_data.process.transcription.rnaData["isTRna"])[0]
		self.idx_rprotein = np.where(sim_data.process.transcription.rnaData['isRProtein'])[0]
		self.idx_rnap = np.where(sim_data.process.transcription.rnaData['isRnap'])[0]

		# Synthesis probabilities for different categories of genes
		self.rnaSynthProbFractions = sim_data.process.transcription.rnaSynthProbFraction
		self.rnaSynthProbRProtein = sim_data.process.transcription.rnaSynthProbRProtein
		self.rnaSynthProbRnaPolymerase = sim_data.process.transcription.rnaSynthProbRnaPolymerase


	def calculateRequest(self):
		# Get all inactive RNA polymerases
		self.inactiveRnaPolys.requestAll()

		# Get attributes of promoters
		promoters = self.promoters.molecules_read_only()
		self.TU_index, bound_TF = promoters.attrs("TU_index", "bound_TF")

		# Read current environment
		current_nutrients = self._external_states['Environment'].nutrients

		if self.full_chromosomes.total_counts()[0] > 0:
			# Calculate probabilities of the RNAP binding to each promoter
			self.promoter_init_probs = (self.basal_prob[self.TU_index] +
				np.multiply(self.delta_prob_matrix[self.TU_index, :], bound_TF).sum(axis=1))

			if len(self.genetic_perturbations) > 0:
				self._rescale_initiation_probs(
					self.genetic_perturbations["fixedRnaIdxs"],
					self.genetic_perturbations["fixedSynthProbs"])

			# Adjust probabilities to not be negative
			self.promoter_init_probs[self.promoter_init_probs < 0] = 0.0
			self.promoter_init_probs /= self.promoter_init_probs.sum()

			# Adjust synthesis probabilities depending on environment
			synthProbFractions = self.rnaSynthProbFractions[current_nutrients]

			# Create masks for different types of RNAs
			is_mrna = np.isin(self.TU_index, self.idx_mrna)
			is_trna = np.isin(self.TU_index, self.idx_trna)
			is_rrna = np.isin(self.TU_index, self.idx_rrna)
			is_rprotein = np.isin(self.TU_index, self.idx_rprotein)
			is_rnap = np.isin(self.TU_index, self.idx_rnap)
			is_fixed = is_trna | is_rrna | is_rprotein | is_rnap

			# Rescale initiation probabilities based on type of RNA
			self.promoter_init_probs[is_mrna] *= synthProbFractions["mRna"] / self.promoter_init_probs[is_mrna].sum()
			self.promoter_init_probs[is_trna] *= synthProbFractions["tRna"] / self.promoter_init_probs[is_trna].sum()
			self.promoter_init_probs[is_rrna] *= synthProbFractions["rRna"] / self.promoter_init_probs[is_rrna].sum()

			# Set fixed synthesis probabilities for RProteins and RNAPs
			self._rescale_initiation_probs(
				np.concatenate((self.idx_rprotein, self.idx_rnap)),
				np.concatenate((
					self.rnaSynthProbRProtein[current_nutrients],
					self.rnaSynthProbRnaPolymerase[current_nutrients]
					)))

			assert self.promoter_init_probs[is_fixed].sum() < 1.0

			# Scale remaining synthesis probabilities accordingly
			scaleTheRestBy = (1. - self.promoter_init_probs[is_fixed].sum()) / self.promoter_init_probs[~is_fixed].sum()
			self.promoter_init_probs[~is_fixed] *= scaleTheRestBy

		# If there are no chromosomes in the cell, set all probs to zero
		else:
			self.promoter_init_probs = np.zeros(len(self.TU_index))

		self.fracActiveRnap = self.fracActiveRnapDict[current_nutrients]
		self.rnaPolymeraseElongationRate = self.rnaPolymeraseElongationRateDict[current_nutrients]


	def evolveState(self):
		# Construct matrix that maps promoters to transcription units
		n_promoters = len(self.TU_index)
		TU_to_promoter = scipy.sparse.csr_matrix(
			(np.ones(n_promoters), (self.TU_index, np.arange(n_promoters))),
			shape = (self.n_TUs, n_promoters))

		# Compute synthesis probabilities of each transcription unit
		TU_synth_probs = TU_to_promoter.dot(self.promoter_init_probs)
		self.writeToListener("RnaSynthProb", "rnaSynthProb", TU_synth_probs)

		# Shuffle synthesis probabilities if we're running the variant that
		# calls this (In general, this should lead to a cell which does not
		# grow and divide)
		if self.shuffleIdxs is not None:
			self._rescale_initiation_probs(
				np.arange(self.n_TUs),
				TU_synth_probs[self.shuffleIdxs])

		# no synthesis if no chromosome
		if self.full_chromosomes.total_counts()[0] == 0:
			return

		# Calculate RNA polymerases to activate based on probabilities
		self.activationProb = self._calculateActivationProb(
			self.fracActiveRnap, self.rnaLengths,
			self.rnaPolymeraseElongationRate, TU_synth_probs)
		rnaPolyToActivate = np.int64(
			self.activationProb * self.inactiveRnaPolys.count())

		if rnaPolyToActivate == 0:
			return

		#### Growth control code ####
		# Sample a multinomial distribution of initiation probabilities to
		# determine what promoters are initialized
		n_initiations = self.randomState.multinomial(rnaPolyToActivate, self.promoter_init_probs)

		# Build list of RNA indexes
		rnaIndexes = np.repeat(self.TU_index, n_initiations)

		# Create the active RNA polymerases
		self.activeRnaPolys.moleculesNew(
			rnaPolyToActivate, rnaIndex = rnaIndexes)
		self.inactiveRnaPolys.countDec(n_initiations.sum())

		# Create masks for ribosomal RNAs
		is_5Srrna = np.isin(self.TU_index, self.idx_5Srrna)
		is_16Srrna = np.isin(self.TU_index, self.idx_16Srrna)
		is_23Srrna = np.isin(self.TU_index, self.idx_23Srrna)

		# Write outputs to listeners
		self.writeToListener(
			"RibosomeData", "rrn16S_produced", n_initiations[is_16Srrna].sum())
		self.writeToListener(
			"RibosomeData", "rrn23S_produced", n_initiations[is_23Srrna].sum())
		self.writeToListener(
			"RibosomeData", "rrn5S_produced", n_initiations[is_5Srrna].sum())

		self.writeToListener(
			"RibosomeData", "rrn16S_init_prob",
			n_initiations[is_16Srrna].sum() / float(rnaPolyToActivate))
		self.writeToListener(
			"RibosomeData", "rrn23S_init_prob",
			n_initiations[is_23Srrna].sum() / float(rnaPolyToActivate))
		self.writeToListener("RibosomeData", "rrn5S_init_prob",
			n_initiations[is_5Srrna].sum() / float(rnaPolyToActivate))

		self.writeToListener("RibosomeData", "total_rna_init", rnaPolyToActivate)

		self.writeToListener("RnapData", "didInitialize", rnaPolyToActivate)
		self.writeToListener("RnapData", "rnaInitEvent", TU_to_promoter.dot(n_initiations))


	def _calculateActivationProb(self, fracActiveRnap, rnaLengths, rnaPolymeraseElongationRate, synthProb):
		"""
		Calculate expected RNAP termination rate based on RNAP elongation rate
		- allTranscriptionTimes: Vector of times required to transcribe each
		transcript
		- allTranscriptionTimestepCounts: Vector of numbers of timesteps
		required to transcribe each transcript
		- averageTranscriptionTimeStepCounts: Average number of timesteps
		required to transcribe a transcript, weighted by synthesis
		probabilities of each transcript
		- expectedTerminationRate: Average number of terminations in one
		timestep for one transcript
		"""
		allTranscriptionTimes = 1. / rnaPolymeraseElongationRate * rnaLengths
		allTranscriptionTimestepCounts = np.ceil(
			(1. / (self.timeStepSec() * units.s) * allTranscriptionTimes).asNumber()
			)
		averageTranscriptionTimestepCounts = np.dot(
			synthProb, allTranscriptionTimestepCounts)
		expectedTerminationRate = 1. / averageTranscriptionTimestepCounts

		"""
		Modify given fraction of active RNAPs to take into account early
		terminations in between timesteps
		- allFractionTimeInactive: Vector of probabilities an "active" RNAP
		will in effect be "inactive" because it has terminated during a
		timestep
		- averageFractionTimeInactive: Average probability of an "active" RNAP
		being in effect "inactive", weighted by synthesis probabilities
		- effectiveFracActiveRnap: New higher "goal" for fraction of active
		RNAP, considering that the "effective" fraction is lower than what the
		listener sees
		"""
		allFractionTimeInactive = 1 - (
			1. / (self.timeStepSec() * units.s) * allTranscriptionTimes).asNumber() / allTranscriptionTimestepCounts
		averageFractionTimeInactive = np.dot(allFractionTimeInactive, synthProb)
		effectiveFracActiveRnap = fracActiveRnap * 1/(1 - averageFractionTimeInactive)

		# Return activation probability that will balance out the expected termination rate
		return effectiveFracActiveRnap * expectedTerminationRate /(1 - effectiveFracActiveRnap)


	def _rescale_initiation_probs(self, fixed_indexes, fixed_synth_probs):
		"""
		Rescales the initiation probabilities of each promoter such that the
		total synthesis probabilities of certain types of RNAs are fixed to
		a predetermined value. For instance, if there are two copies of
		promoters for RNA A, whose synthesis probability should be fixed to
		0.1, each promoter is given an initiation probability of 0.05.
		"""
		for idx, synth_prob in izip(fixed_indexes, fixed_synth_probs):
			fixed_mask = (self.TU_index == idx)
			self.promoter_init_probs[fixed_mask] = synth_prob / fixed_mask.sum()

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

import itertools

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

		# Create recruitment matrix for transcription regulation
		recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames
		recruitmentData = sim_data.process.transcription_regulation.recruitmentData
		self.recruitmentMatrix = scipy.sparse.csr_matrix(
			(recruitmentData["hV"], (recruitmentData["hI"], recruitmentData["hJ"])),
			shape = recruitmentData["shape"]
			)

		self.maxRibosomeElongationRate = float(
			sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s))

		# Determine changes from genetic perturbations
		self.genetic_perturbations = {}
		perturbations = getattr(sim_data, "genetic_perturbations", {})

		if len(perturbations) > 0:
			probability_indexes = [
				(index, sim_data.genetic_perturbations[rna_data['id']])
				for index, rna_data in
				enumerate(sim_data.process.transcription.rnaData)
				if rna_data['id'] in sim_data.genetic_perturbations]

			self.genetic_perturbations = {
				'fixedRnaIdxs': map(lambda pair: pair[0], probability_indexes),
				'fixedSynthProbs': map(lambda pair: pair[1], probability_indexes)
				}

		# If initiationShuffleIdxs does not exist, set value to None
		self.shuffleIdxs = getattr(
			sim_data.process.transcription, 'initiationShuffleIdxs', None)

		# Views
		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')
		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")
		self.full_chromosomes = self.uniqueMoleculesView('fullChromosome')
		self.recruitmentView = self.bulkMoleculesView(recruitmentColNames)

		# ID Groups
		self.is_16SrRNA = sim_data.process.transcription.rnaData['isRRna16S']
		self.is_23SrRNA = sim_data.process.transcription.rnaData['isRRna23S']
		self.is_5SrRNA = sim_data.process.transcription.rnaData['isRRna5S']
		self.isRRna = sim_data.process.transcription.rnaData['isRRna']
		self.isMRna = sim_data.process.transcription.rnaData["isMRna"]
		self.isTRna = sim_data.process.transcription.rnaData["isTRna"]
		self.isRProtein = sim_data.process.transcription.rnaData['isRProtein']
		self.isRnap = sim_data.process.transcription.rnaData['isRnap']
		self.setIdxs = self.isRRna | self.isTRna | self.isRProtein | self.isRnap

		# Synthesis probabilities for different categories of genes
		self.rnaSynthProbFractions = sim_data.process.transcription.rnaSynthProbFraction
		self.rnaSynthProbRProtein = sim_data.process.transcription.rnaSynthProbRProtein
		self.rnaSynthProbRnaPolymerase = sim_data.process.transcription.rnaSynthProbRnaPolymerase


	def calculateRequest(self):
		# Get all inactive RNA polymerases
		self.inactiveRnaPolys.requestAll()

		# Read current environment
		current_nutrients = self._external_states['Environment'].nutrients

		if self.full_chromosomes.total_counts()[0] > 0:
			# Calculate synthesis probabilities based on transcription regulation
			self.rnaSynthProb = self.recruitmentMatrix.dot(
				self.recruitmentView.total_counts())
			if len(self.genetic_perturbations) > 0:
				self.rnaSynthProb[self.genetic_perturbations["fixedRnaIdxs"]] = self.genetic_perturbations["fixedSynthProbs"]

			# Adjust probabilities to not be negative
			self.rnaSynthProb[self.rnaSynthProb < 0] = 0.0
			self.rnaSynthProb /= self.rnaSynthProb.sum()

			# Adjust synthesis probabilities depending on environment
			synthProbFractions = self.rnaSynthProbFractions[current_nutrients]

			# Allocate synthesis probabilities based on type of RNA
			self.rnaSynthProb[self.isMRna] *= synthProbFractions["mRna"] / self.rnaSynthProb[self.isMRna].sum()
			self.rnaSynthProb[self.isTRna] *= synthProbFractions["tRna"] / self.rnaSynthProb[self.isTRna].sum()
			self.rnaSynthProb[self.isRRna] *= synthProbFractions["rRna"] / self.rnaSynthProb[self.isRRna].sum()

			# Set fixed synthesis probabilities for RProteins and RNAPs
			self.rnaSynthProb[self.isRProtein] = self.rnaSynthProbRProtein[current_nutrients]
			self.rnaSynthProb[self.isRnap] = self.rnaSynthProbRnaPolymerase[current_nutrients]

			assert self.rnaSynthProb[self.setIdxs].sum() < 1.0

			# Scale remaining synthesis probabilities accordingly
			scaleTheRestBy = (1. - self.rnaSynthProb[self.setIdxs].sum()) / self.rnaSynthProb[~self.setIdxs].sum()
			self.rnaSynthProb[~self.setIdxs] *= scaleTheRestBy

			# Shuffle initiation rates if we're running the variant that calls this
			# (In general, this should lead to a cell which does not grow and
			# divide)
			if self.shuffleIdxs is not None:
				self.rnaSynthProb = self.rnaSynthProb[self.shuffleIdxs]

		# If there are no chromosomes in the cell, set all probs to zero
		else:
			self.rnaSynthProb = np.zeros(self.recruitmentMatrix.shape[0])

		self.fracActiveRnap = self.fracActiveRnapDict[current_nutrients]
		self.rnaPolymeraseElongationRate = self.rnaPolymeraseElongationRateDict[current_nutrients]


	def evolveState(self):
		self.writeToListener("RnaSynthProb", "rnaSynthProb", self.rnaSynthProb)

		# no synthesis if no chromosome
		if self.full_chromosomes.total_counts()[0] == 0:
			return

		# Calculate RNA polymerases to activate based on probabilities
		self.activationProb = self._calculateActivationProb(
			self.fracActiveRnap, self.rnaLengths,
			self.rnaPolymeraseElongationRate, self.rnaSynthProb)
		rnaPolyToActivate = np.int64(
			self.activationProb * self.inactiveRnaPolys.count())

		if rnaPolyToActivate == 0:
			return

		#### Growth control code ####
		# Sample a multinomial distribution of synthesis probabilities to
		# determine what molecules are initialized
		nNewRnas = self.randomState.multinomial(rnaPolyToActivate, self.rnaSynthProb)

		# Build list of RNA indexes
		rnaIndexes = np.empty(rnaPolyToActivate, np.int64)
		startIndex = 0
		nonzeroCount = (nNewRnas > 0)
		for rnaIndex, counts in itertools.izip(
				np.arange(nNewRnas.size)[nonzeroCount], nNewRnas[nonzeroCount]):
			rnaIndexes[startIndex:startIndex+counts] = rnaIndex
			startIndex += counts

		# Create the active RNA polymerases
		activeRnaPolys = self.activeRnaPolys.moleculesNew(
			"activeRnaPoly", rnaPolyToActivate)
		activeRnaPolys.attrIs(rnaIndex = rnaIndexes)
		self.inactiveRnaPolys.countDec(nNewRnas.sum())

		# Write outputs to listeners
		self.writeToListener(
			"RibosomeData", "rrn16S_produced", nNewRnas[self.is_16SrRNA].sum())
		self.writeToListener(
			"RibosomeData", "rrn23S_produced", nNewRnas[self.is_23SrRNA].sum())
		self.writeToListener(
			"RibosomeData", "rrn5S_produced", nNewRnas[self.is_5SrRNA].sum())

		self.writeToListener(
			"RibosomeData", "rrn16S_init_prob",
			nNewRnas[self.is_16SrRNA].sum() / float(nNewRnas.sum())
			)
		self.writeToListener(
			"RibosomeData", "rrn23S_init_prob",
			nNewRnas[self.is_23SrRNA].sum() / float(nNewRnas.sum())
			)
		self.writeToListener("RibosomeData", "rrn5S_init_prob",
			nNewRnas[self.is_5SrRNA].sum() / float(nNewRnas.sum())
			)

		self.writeToListener("RibosomeData", "total_rna_init", nNewRnas.sum())

		self.writeToListener("RnapData", "didInitialize", nNewRnas.sum())
		self.writeToListener("RnapData", "rnaInitEvent", nNewRnas)

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

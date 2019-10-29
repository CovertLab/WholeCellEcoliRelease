#!/usr/bin/env python

"""
TranscriptInitiation

Transcription initiation sub-model.

TODO:
- use transcription units instead of single genes
- match sigma factors to promoters
- implement transcriptional regulation
- modulate initiation probabilities as a function of gene copy number
- match measured levels of active RNA polymerase instead of initiating to completion

@author: John Mason
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

		# Determine changes from genetic perturbations
		self.genetic_perturbations = {}
		perturbations = getattr(sim_data, "genetic_perturbations", {})
		if len(perturbations) > 0:
			rnaIdxs, synthProbs = zip(*[(int(np.where(sim_data.process.transcription.rnaData["id"] == rnaId)[0]), synthProb) for rnaId, synthProb in sim_data.genetic_perturbations.iteritems()])
			fixedSynthProbs = [synthProb for (rnaIdx, synthProb) in sorted(zip(rnaIdxs, synthProbs), key = lambda pair: pair[0])]
			fixedRnaIdxs = [rnaIdx for (rnaIdx, synthProb) in sorted(zip(rnaIdxs, synthProbs), key = lambda pair: pair[0])]
			self.genetic_perturbations = {"fixedRnaIdxs": fixedRnaIdxs, "fixedSynthProbs": fixedSynthProbs}

		# If initiationShuffleIdxs does not exist, set value to None
		self.shuffleIdxs = getattr(sim_data.process.transcription, 'initiationShuffleIdxs', None)

		# Views
		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')
		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")
		self.chromosomes = self.bulkMoleculeView('CHROM_FULL[c]')
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
		self.isRegulated = np.array([1 if x[:-3] in sim_data.process.transcription_regulation.targetTf or x in perturbations else 0 for x in sim_data.process.transcription.rnaData["id"]], dtype = np.bool)
		self.setIdxs = self.isRRna | self.isTRna | self.isRProtein | self.isRnap | self.isRegulated

		# Synthesis probabilities for different categories of genes
		self.rnaSynthProbFractions = sim_data.process.transcription.rnaSynthProbFraction
		self.rnaSynthProbRProtein = sim_data.process.transcription.rnaSynthProbRProtein
		self.rnaSynthProbRnaPolymerase = sim_data.process.transcription.rnaSynthProbRnaPolymerase

	def calculateRequest(self):
		# Get all inactive RNA polymerases
		self.inactiveRnaPolys.requestAll()

		# Calculate synthesis probabilities based on transcription regulation
		self.rnaSynthProb = self.recruitmentMatrix.dot(self.recruitmentView.total())
		if len(self.genetic_perturbations) > 0:
			self.rnaSynthProb[self.genetic_perturbations["fixedRnaIdxs"]] = self.genetic_perturbations["fixedSynthProbs"]
		regProbs = self.rnaSynthProb[self.isRegulated]

		# Adjust probabilities to not be negative
		self.rnaSynthProb[self.rnaSynthProb < 0] = 0.0
		self.rnaSynthProb /= self.rnaSynthProb.sum()
		if np.any(self.rnaSynthProb < 0):
			raise Exception("Have negative RNA synthesis probabilities")

		# Adjust synthesis probabilities depending on environment
		current_nutrients = self._external_states['Environment'].nutrients

		synthProbFractions = self.rnaSynthProbFractions[current_nutrients]
		self.rnaSynthProb[self.isMRna] *= synthProbFractions["mRna"] / self.rnaSynthProb[self.isMRna].sum()
		self.rnaSynthProb[self.isTRna] *= synthProbFractions["tRna"] / self.rnaSynthProb[self.isTRna].sum()
		self.rnaSynthProb[self.isRRna] *= synthProbFractions["rRna"] / self.rnaSynthProb[self.isRRna].sum()
		self.rnaSynthProb[self.isRegulated] = regProbs
		self.rnaSynthProb[self.isRProtein] = self.rnaSynthProbRProtein[current_nutrients]
		self.rnaSynthProb[self.isRnap] = self.rnaSynthProbRnaPolymerase[current_nutrients]
		self.rnaSynthProb[self.rnaSynthProb < 0] = 0
		scaleTheRestBy = (1. - self.rnaSynthProb[self.setIdxs].sum()) / self.rnaSynthProb[~self.setIdxs].sum()
		self.rnaSynthProb[~self.setIdxs] *= scaleTheRestBy

		# Shuffle initiation rates if we're running the variant that calls this
		# (In general, this should lead to a cell which does not grow and divide)
		if self.shuffleIdxs is not None:
			self.rnaSynthProb = self.rnaSynthProb[self.shuffleIdxs]

		self.fracActiveRnap = self.fracActiveRnapDict[current_nutrients]
		self.rnaPolymeraseElongationRate = self.rnaPolymeraseElongationRateDict[current_nutrients]

	def evolveState(self):
		self.writeToListener("RnaSynthProb", "rnaSynthProb", self.rnaSynthProb)

		# no synthesis if no chromosome
		if self.chromosomes.total()[0] == 0:
			return

		# Calculate RNA polymerases to activate based on probabilities
		self.activationProb = self._calculateActivationProb(self.fracActiveRnap, self.rnaLengths, self.rnaPolymeraseElongationRate, self.rnaSynthProb)
		rnaPolyToActivate = np.int64(self.activationProb * self.inactiveRnaPolys.count())
		if rnaPolyToActivate == 0:
			return

		#### Growth control code ####
		# Sample a multinomial distribution of synthesis probabilities to determine what molecules are initialized
		nNewRnas = self.randomState.multinomial(rnaPolyToActivate, self.rnaSynthProb)

		# Build list of RNA indexes
		rnaIndexes = np.empty(rnaPolyToActivate, np.int64)
		startIndex = 0
		nonzeroCount = (nNewRnas > 0)
		for rnaIndex, counts in itertools.izip(np.arange(nNewRnas.size)[nonzeroCount], nNewRnas[nonzeroCount]):
			rnaIndexes[startIndex:startIndex+counts] = rnaIndex
			startIndex += counts

		# Create the active RNA polymerases
		activeRnaPolys = self.activeRnaPolys.moleculesNew("activeRnaPoly", rnaPolyToActivate)
		activeRnaPolys.attrIs(rnaIndex = rnaIndexes)
		self.inactiveRnaPolys.countDec(nNewRnas.sum())

		# Write outputs to listeners
		self.writeToListener("RibosomeData", "rrn16S_produced", nNewRnas[self.is_16SrRNA].sum())
		self.writeToListener("RibosomeData", "rrn23S_produced", nNewRnas[self.is_23SrRNA].sum())
		self.writeToListener("RibosomeData", "rrn5S_produced", nNewRnas[self.is_5SrRNA].sum())

		self.writeToListener("RibosomeData", "rrn16S_init_prob", nNewRnas[self.is_16SrRNA].sum() / float(nNewRnas.sum()))
		self.writeToListener("RibosomeData", "rrn23S_init_prob", nNewRnas[self.is_23SrRNA].sum() / float(nNewRnas.sum()))
		self.writeToListener("RibosomeData", "rrn5S_init_prob", nNewRnas[self.is_5SrRNA].sum() / float(nNewRnas.sum()))

		self.writeToListener("RibosomeData", "total_rna_init", nNewRnas.sum())

		self.writeToListener("RnapData", "didInitialize", nNewRnas.sum())
		self.writeToListener("RnapData", "rnaInitEvent", nNewRnas)

	def _calculateActivationProb(self, fracActiveRnap, rnaLengths, rnaPolymeraseElongationRate, synthProb):
		'''
		Edge case: 100% RNA polymerase activity (relevant to paper investigations).
		Return: 1
		'''
		if fracActiveRnap == 1:
			return 1.

		''' Calculate expected RNAP termination rate based on RNAP elongation rate
		allTranscriptionTimes: Vector of times required to transcribe each transcript
		allTranscriptionTimestepCounts: Vector of numbers of timesteps required to transcribe each transcript
		averageTranscriptionTimeStepCounts: Average number of timesteps required to transcribe a transcript, weighted by synthesis probabilities
		expectedTerminationRate: Average number of terminations in one timestep for one transcript
		'''
		allTranscriptionTimes = 1. / rnaPolymeraseElongationRate * rnaLengths
		allTranscriptionTimestepCounts = np.ceil((1. / (self.timeStepSec() * units.s) * allTranscriptionTimes).asNumber())
		averageTranscriptionTimestepCounts = np.dot(synthProb, allTranscriptionTimestepCounts)
		expectedTerminationRate = 1. / averageTranscriptionTimestepCounts

		''' Modify given fraction of active RNAPs to take into account early terminations in between timesteps
		allFractionTimeInactive: Vector of probabilities an "active" RNAP will in effect be "inactive" because it has terminated during a timestep
		averageFractionTimeInactive: Average probability of an "active" RNAP being in effect "inactive", weighted by synthesis probabilities
		effectiveFracActiveRnap: New higher "goal" for fraction of active RNAP, considering that the "effective" fraction is lower than what the listener sees
		'''
		allFractionTimeInactive = 1 - ( 1. / (self.timeStepSec() * units.s) * allTranscriptionTimes).asNumber() / allTranscriptionTimestepCounts
		averageFractionTimeInactive = np.dot(allFractionTimeInactive, synthProb)
		effectiveFracActiveRnap = fracActiveRnap * 1 / (1 - averageFractionTimeInactive)

		# Return activation probability that will balance out the expected termination rate
		return effectiveFracActiveRnap * expectedTerminationRate / (1 - effectiveFracActiveRnap)

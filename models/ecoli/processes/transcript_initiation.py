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


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(TranscriptInitiation, self).initialize(sim, sim_data)

		# Load parameters

		self.fracActiveRnapDict = sim_data.process.transcription.rnapFractionActiveDict

		self.rnaLengths = sim_data.process.transcription.rnaData["length"]

		self.rnaPolymeraseElongationRateDict = sim_data.process.transcription.rnaPolymeraseElongationRateDict

		self.rnaSynthProb = None

		recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames

		recruitmentData = sim_data.process.transcription_regulation.recruitmentData
		self.recruitmentMatrix = scipy.sparse.csr_matrix(
				(recruitmentData["hV"], (recruitmentData["hI"], recruitmentData["hJ"])),
				shape = recruitmentData["shape"]
			)
		self.tfsBound = None

		self.maxRibosomeElongationRate = float(sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s))
		self.is_16SrRNA = sim_data.process.transcription.rnaData['isRRna16S']
		self.is_23SrRNA = sim_data.process.transcription.rnaData['isRRna23S']
		self.is_5SrRNA = sim_data.process.transcription.rnaData['isRRna5S']
		self.is_mrRNA = sim_data.process.transcription.rnaData['isMRna']

		self.genetic_perturbations = {}
		perturbations = {}
		if hasattr(sim_data, "genetic_perturbations") and sim_data.genetic_perturbations != None and len(sim_data.genetic_perturbations) > 0:
			rnaIdxs, synthProbs = zip(*[(int(np.where(sim_data.process.transcription.rnaData["id"] == rnaId)[0]), synthProb) for rnaId, synthProb in sim_data.genetic_perturbations.iteritems()])
			fixedSynthProbs = [synthProb for (rnaIdx, syntheProb) in sorted(zip(rnaIdxs, synthProbs), key = lambda pair: pair[0])]
			fixedRnaIdxs = [rnaIdx for (rnaIdx, syntheProb) in sorted(zip(rnaIdxs, synthProbs), key = lambda pair: pair[0])]
			self.genetic_perturbations = {"fixedRnaIdxs": fixedRnaIdxs, "fixedSynthProbs": fixedSynthProbs}
			perturbations = sim_data.genetic_perturbations

		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')

		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")

		self.chromosomes = self.bulkMoleculeView('CHROM_FULL[c]')
		
		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')

		self.r_protein = self.bulkMoleculesView(sim_data.moleculeGroups.rProteins)

		# ID Groups

		self.is_16SrRNA = sim_data.process.transcription.rnaData['isRRna16S']
		self.is_23SrRNA = sim_data.process.transcription.rnaData['isRRna23S']
		self.is_5SrRNA = sim_data.process.transcription.rnaData['isRRna5S']

		self.recruitmentView = self.bulkMoleculesView(recruitmentColNames)

		self.isRRna = sim_data.process.transcription.rnaData['isRRna']
		self.isMRna = sim_data.process.transcription.rnaData["isMRna"]
		self.isTRna = sim_data.process.transcription.rnaData["isTRna"]
		self.isRProtein = sim_data.process.transcription.rnaData['isRProtein']
		self.isRnap = sim_data.process.transcription.rnaData['isRnap']
		self.notPolymerase = np.logical_and(np.logical_and(np.logical_not(self.isRRna),np.logical_not(self.isRProtein)), np.logical_not(self.isRnap))
		self.isRegulated = np.array([1 if x[:-3] in sim_data.process.transcription_regulation.targetTf or x in perturbations else 0 for x in sim_data.process.transcription.rnaData["id"]], dtype = np.bool)
		self.setIdxs = self.isRRna | self.isTRna | self.isRProtein | self.isRnap | self.isRegulated

		assert (self.isRRna + self.isRProtein + self.isRnap + self.notPolymerase).sum() == self.rnaLengths.asNumber().size

		self.rProteinToRRnaRatioVector = None
		self.rnaSynthProbFractions = sim_data.process.transcription.rnaSynthProbFraction
		self.rnaSynthProbRProtein = sim_data.process.transcription.rnaSynthProbRProtein
		self.rnaSynthProbRnaPolymerase = sim_data.process.transcription.rnaSynthProbRnaPolymerase

	def calculateRequest(self):
		self.inactiveRnaPolys.requestAll()
		self.rnaSynthProb = self.recruitmentMatrix.dot(self.recruitmentView.total())
		if len(self.genetic_perturbations) > 0:
			self.rnaSynthProb[self.genetic_perturbations["fixedRnaIdxs"]] = self.genetic_perturbations["fixedSynthProbs"]
		regProbs = self.rnaSynthProb[self.isRegulated]
		self.rnaSynthProb[self.rnaSynthProb < 0] = 0.
		self.rnaSynthProb /= self.rnaSynthProb.sum()
		if np.any(self.rnaSynthProb < 0):
			raise Exception, "Have negative RNA synthesis probabilities"

		assert np.allclose(self.rnaSynthProb.sum(),1.)
		assert np.all(self.rnaSynthProb >= 0.)

		synthProbFractions = self.rnaSynthProbFractions[self._sim.processes["PolypeptideElongation"].currentNutrients]
		self.rnaSynthProb[self.isMRna] *= synthProbFractions["mRna"] / self.rnaSynthProb[self.isMRna].sum()
		self.rnaSynthProb[self.isTRna] *= synthProbFractions["tRna"] / self.rnaSynthProb[self.isTRna].sum()
		self.rnaSynthProb[self.isRRna] *= synthProbFractions["rRna"] / self.rnaSynthProb[self.isRRna].sum()
		self.rnaSynthProb[self.isRegulated] = regProbs
		self.rnaSynthProb[self.isRProtein] = self.rnaSynthProbRProtein[self._sim.processes["PolypeptideElongation"].currentNutrients]
		self.rnaSynthProb[self.isRnap] = self.rnaSynthProbRnaPolymerase[self._sim.processes["PolypeptideElongation"].currentNutrients]
		self.rnaSynthProb[self.rnaSynthProb < 0] = 0
		scaleTheRestBy = (1. - self.rnaSynthProb[self.setIdxs].sum()) / self.rnaSynthProb[~self.setIdxs].sum()
		self.rnaSynthProb[~self.setIdxs] *= scaleTheRestBy

		assert np.allclose(self.rnaSynthProb.sum(),1.)
		assert np.all(self.rnaSynthProb >= 0.)

		self.fracActiveRnap = self.fracActiveRnapDict[self._sim.processes["PolypeptideElongation"].currentNutrients]
		self.rnaPolymeraseElongationRate = self.rnaPolymeraseElongationRateDict[self._sim.processes["PolypeptideElongation"].currentNutrients]

		# self.rProteinToRRnaRatioVector = self.rnaSynthProb[self.isRProtein] / self.rnaSynthProb[self.isRRna][0]

	# Calculate temporal evolution
	def evolveState(self):

		self.writeToListener("RnaSynthProb", "rnaSynthProb", self.rnaSynthProb)

		# no synthesis if no chromosome
		if self.chromosomes.total()[0] == 0:
			return

		self.activationProb = self._calculateActivationProb(
			self.fracActiveRnap,
			self.rnaLengths,
			self.rnaPolymeraseElongationRate,
			self.rnaSynthProb,
			)

		# Sample a multinomial distribution of synthesis probabilities to 
		# determine what molecules are initialized

		inactiveRnaPolyCount = self.inactiveRnaPolys.count()

		rnaPolyToActivate = np.int64(self.activationProb * inactiveRnaPolyCount)

		if rnaPolyToActivate == 0:
			return

		#### Growth control code ####
		nNewRnas = self.randomState.multinomial(rnaPolyToActivate,
			self.rnaSynthProb)

		self.writeToListener("RibosomeData", "rrn16S_produced", nNewRnas[self.is_16SrRNA].sum())
		self.writeToListener("RibosomeData", "rrn23S_produced", nNewRnas[self.is_23SrRNA].sum())		
		self.writeToListener("RibosomeData", "rrn5S_produced", nNewRnas[self.is_5SrRNA].sum())

		self.writeToListener("RibosomeData", "rrn16S_init_prob", nNewRnas[self.is_16SrRNA].sum() / float(nNewRnas.sum()))
		self.writeToListener("RibosomeData", "rrn23S_init_prob", nNewRnas[self.is_23SrRNA].sum() / float(nNewRnas.sum()))
		self.writeToListener("RibosomeData", "rrn5S_init_prob", nNewRnas[self.is_5SrRNA].sum() / float(nNewRnas.sum()))

		self.writeToListener("RibosomeData", "total_rna_init", nNewRnas.sum())

		nonzeroCount = (nNewRnas > 0)

		assert nNewRnas.sum() == rnaPolyToActivate

		# Build list of RNA indexes

		rnaIndexes = np.empty(rnaPolyToActivate, np.int64)

		startIndex = 0
		for rnaIndex, counts in itertools.izip(
				np.arange(nNewRnas.size)[nonzeroCount],
				nNewRnas[nonzeroCount]
				):

			rnaIndexes[startIndex:startIndex+counts] = rnaIndex

			startIndex += counts

		# Create the active RNA polymerases

		activeRnaPolys = self.activeRnaPolys.moleculesNew(
			"activeRnaPoly",
			rnaPolyToActivate
			)

		activeRnaPolys.attrIs(
			rnaIndex = rnaIndexes
			)

		self.inactiveRnaPolys.countDec(nNewRnas.sum())

		self.writeToListener("RnapData", "didInitialize", nNewRnas.sum())
		self.writeToListener("RnapData", "rnaInitEvent", nNewRnas)

	def _calculateActivationProb(self, fracActiveRnap, rnaLengths, rnaPolymeraseElongationRate, synthProb):
		expectedTranscriptionTime = 1. / rnaPolymeraseElongationRate * rnaLengths

		expectedTranscriptionTimesteps = np.ceil(
			(1. / (self.timeStepSec() * units.s) * expectedTranscriptionTime).asNumber()
			)

		averageTranscriptionTimesteps = np.dot(synthProb, expectedTranscriptionTimesteps)

		expectedTerminationRate = 1. / averageTranscriptionTimesteps

		expectedFractionTimeInactive = np.dot(
			1 - ( 1. / (self.timeStepSec() * units.s) * expectedTranscriptionTime).asNumber() / expectedTranscriptionTimesteps,
			synthProb
			)

		effectiveFractionActive = fracActiveRnap * 1 / (1 - expectedFractionTimeInactive)

		return effectiveFractionActive * expectedTerminationRate / (1 - effectiveFractionActive)

	def calculateRrnInitRate(self, rrn_count, elngRate):
		'''
		Returns total initiation rate of rRNA across all promoters
		In units of initiations / s / fg
		'''
		fitInitiationRate = 0.0168 * np.exp(-0.272 * (self.maxRibosomeElongationRate - elngRate))

		return (1 / units.s / units.fg) * fitInitiationRate

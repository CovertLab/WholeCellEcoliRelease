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
		# Parameters
		self.rnaSynthProb = None

		# Views
		self.activeRnaPolys = None
		self.inactiveRnaPolys = None

		super(TranscriptInitiation, self).__init__()


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(TranscriptInitiation, self).initialize(sim, sim_data)

		# Load parameters

		self.fracActiveRnap = sim_data.growthRateParameters.fractionActiveRnap

		self.rnaLengths = sim_data.process.transcription.rnaData["length"]

		self.rnaPolymeraseElongationRate = sim_data.growthRateParameters.rnaPolymeraseElongationRate

		self.rnaSynthProb = None

		recruitmentColNames = sim_data.process.transcription_regulation.recruitmentColNames

		recruitmentData = sim_data.process.transcription_regulation.recruitmentData
		self.recruitmentMatrix = scipy.sparse.csr_matrix(
				(recruitmentData["hV"], (recruitmentData["hI"], recruitmentData["hJ"])),
				shape = recruitmentData["shape"]
			)
		self.tfsBound = None

		self.is_16SrRNA = sim_data.process.transcription.rnaData['isRRna16S']
		self.is_23SrRNA = sim_data.process.transcription.rnaData['isRRna23S']
		self.is_5SrRNA = sim_data.process.transcription.rnaData['isRRna5S']


		import copy
		self.rnaSynthProbStandard = copy.copy(self.rnaSynthProb)

		self.maxRibosomeElongationRate = float(sim_data.constants.ribosomeElongationRateMax.asNumber(units.aa / units.s))

		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')

		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")

		self.chromosomes = self.bulkMoleculeView('CHROM_FULL[c]')
		
		self.rrn_operon = self.bulkMoleculeView("rrn_operon")

		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')

		# ID Groups

		self.is_16SrRNA = sim_data.process.transcription.rnaData['isRRna16S']
		self.is_23SrRNA = sim_data.process.transcription.rnaData['isRRna23S']
		self.is_5SrRNA = sim_data.process.transcription.rnaData['isRRna5S']

		self.recruitmentView = self.bulkMoleculesView(recruitmentColNames)

		self.isRRna = sim_data.process.transcription.rnaData['isRRna']
		self.isRProtein = sim_data.process.transcription.rnaData['isRProtein']
		self.isRnap = sim_data.process.transcription.rnaData['isRnap']
		self.notPolymerase = np.logical_and(np.logical_and(np.logical_not(self.isRRna),np.logical_not(self.isRProtein)), np.logical_not(self.isRnap))

		assert (self.isRRna + self.isRProtein + self.isRnap + self.notPolymerase).sum() == self.rnaSynthProb.size

		self.rProteinToRRnaRatioVector = self.rnaSynthProbStandard[self.isRProtein] / self.rnaSynthProbStandard[self.isRRna][0]

		## VARIANT CODE ##
		# self.scaling_factor = sim_data.scaling_factor
		self.scaling_factor = 10
		## VARIANT CODE ##

	def calculateRequest(self):
		self.inactiveRnaPolys.requestAll()
		self.rnaSynthProb = self.recruitmentMatrix.dot(self.recruitmentView.total())
		self.rnaSynthProb /= self.rnaSynthProb.sum()

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

		ribosomeElongationRate = self.readFromListener("RibosomeData", "effectiveElongationRate")
		expectedRibosomeInitiationRate = self.calculateRrnInitRate(self.rrn_operon.total(), ribosomeElongationRate)
		rRnaSynthesisProb = expectedRibosomeInitiationRate.asNumber(1/units.s) * self.timeStepSec() / rnaPolyToActivate
		rProteinSynthesisProb = self.rProteinToRRnaRatioVector * rRnaSynthesisProb

		totalRnapCount = self.activeRnaPolys.total() + self.inactiveRnaPolys.total() or np.array([1])
		totalRibosomeCount = self.activeRibosomes.total() or np.array([1])

		ratioRNAPToRibosome = totalRnapCount / totalRibosomeCount.astype(np.float)
		offset = np.clip(0.25 - ratioRNAPToRibosome, -1 * self.rnaSynthProbStandard[self.isRnap].min(), 1.)
		rnapSynthProb = self.rnaSynthProbStandard[self.isRnap] + (offset / 10)

		print "Expected init rate: {}".format(expectedRibosomeInitiationRate.asNumber(1/units.s))
		self.writeToListener("RibosomeData", "expectedInitRate", expectedRibosomeInitiationRate.asNumber(1/units.s))

		totalRRnaSynthProb = (np.ceil(self.rnaSynthProb[self.isRRna]).sum() * rRnaSynthesisProb) # HACK: Only getting used ribosome rrn operons using this ceil function need to fix this.
		totalRProteinSynthProb = rProteinSynthesisProb.sum()
		totalRnapSynthProb = rnapSynthProb.sum()

		totalPolymeraseComponent = totalRRnaSynthProb + totalRProteinSynthProb + totalRnapSynthProb


		while totalPolymeraseComponent > 1.:
			rRnaSynthesisProb = rRnaSynthesisProb / totalPolymeraseComponent
			rProteinSynthesisProb = rProteinSynthesisProb / totalPolymeraseComponent
			rnapSynthProb = rnapSynthProb / totalPolymeraseComponent

			totalRRnaSynthProb = totalRRnaSynthProb / totalPolymeraseComponent
			totalRProteinSynthProb = totalRProteinSynthProb / totalPolymeraseComponent
			totalRnapSynthProb = totalRnapSynthProb / totalPolymeraseComponent

			totalPolymeraseComponent = totalRRnaSynthProb + totalRProteinSynthProb + totalRnapSynthProb

		print "Correct prob: {}".format(self.rnaSynthProbStandard[self.isRRna].sum() + self.rnaSynthProbStandard[self.isRProtein].sum())
		print "Actual prob: {}".format(totalPolymeraseComponent)

		self.rnaSynthProb[self.isRRna] = rRnaSynthesisProb * np.ceil(self.rnaSynthProb[self.isRRna]) # HACK ALERT: Only getting used ribosome rrn operons using this ceil function need to fix this.

		self.rnaSynthProb[self.isRProtein] = rProteinSynthesisProb

		self.rnaSynthProb[self.isRnap] = rnapSynthProb

		self.rnaSynthProb[self.notPolymerase] = (1 - totalPolymeraseComponent) / self.rnaSynthProbStandard[self.notPolymerase].sum() * self.rnaSynthProbStandard[self.notPolymerase]

		if not np.allclose(self.rnaSynthProb.sum(),1.):
			import ipdb; ipdb.set_trace()
		if not np.all(self.rnaSynthProb >= 0.):
			import ipdb; ipdb.set_trace()
		assert np.allclose(self.rnaSynthProb.sum(),1.)
		assert np.all(self.rnaSynthProb >= 0.)

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
		In units of initiations / min
		'''
		fitInitiationRate = rrn_count[0] * 151.595 * np.exp(0.038*-0.298 * (self.maxRibosomeElongationRate - elngRate)) / self.scaling_factor
		# fitInitiationRate = 151.595 * np.exp(0.038*-0.298 * (self.maxRibosomeElongationRate - elngRate)) / 10.

	#	fitInitiationRate = rrn_count * 151.595 * np.exp(0.038*-0.298 * (self.maxRibosomeElongationRate - elngRate))

		return (1 / units.min) * fitInitiationRate

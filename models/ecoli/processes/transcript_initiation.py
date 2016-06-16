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


		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')

		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")

		self.chromosomes = self.bulkMoleculeView('CHROM_FULL[c]')

		self.recruitmentView = self.bulkMoleculesView(recruitmentColNames)


	def calculateRequest(self):
		self.inactiveRnaPolys.requestAll()
		self.rnaSynthProb = self.recruitmentMatrix.dot(self.recruitmentView.total())
		self.rnaSynthProb /= self.rnaSynthProb.sum()


	# Calculate temporal evolution
	def evolveState(self):

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
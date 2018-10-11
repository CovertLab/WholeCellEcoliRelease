#!/usr/bin/env python

"""
RibosomeData

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/21/14
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

# from numpy.lib.recfunctions import merge_arrays

VERBOSE = False

class RibosomeData(wholecell.listeners.listener.Listener):
	""" RibosomeData """

	_name = 'RibosomeData'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RibosomeData, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RibosomeData, self).initialize(sim, sim_data)

		self.nMonomers = len(sim_data.process.translation.monomerData)

		# Logged quantities
		self.registerLoggedQuantity(
			"Fraction\nribosomes\nstalled",
			"fractionStalled",
			".3f"
			)


	# Allocate memory
	def allocate(self):
		super(RibosomeData, self).allocate()

		# Computed, saved attributes
		self.stallingRateTotal = 0
		self.stallingRateMean = 0
		self.stallingRateStd = 0
		self.fractionStalled = 0

		# Attributes broadcast by the PolypeptideElongation process
		self.ribosomeStalls = np.zeros(0, np.int64)
		self.aaCountInSequence = np.zeros(21, np.int64)
		self.aaCounts = np.zeros(21, np.int64)
		self.actualElongations = 0
		self.actualElongationHist = np.zeros(22, np.int64)
		self.elongationsNonTerminatingHist = np.zeros(22, np.int64)
		self.expectedElongations = 0
		self.didTerminate = 0
		self.didInitialize = 0
		self.terminationLoss = 0
		self.effectiveElongationRate = 0.
		self.rrn16S_produced = 0
		self.rrn23S_produced = 0
		self.rrn5S_produced = 0
		self.rrn16S_init_prob = 0.
		self.rrn23S_init_prob = 0.
		self.rrn5S_init_prob = 0.
		self.total_rna_init = 0
		self.processElongationRate = 0.
		self.translationSupply = np.zeros(21, np.float64)
		self.numTrpATerminated = 0.
		self.probTranslationPerTranscript = np.zeros(self.nMonomers, np.float64)


	def update(self):
		if self.ribosomeStalls.size:
			# TODO: divide rates by time step length
			self.stallingRateTotal = self.ribosomeStalls.sum()
			self.stallingRateMean = self.ribosomeStalls.mean()
			self.stallingRateStd = self.ribosomeStalls.std()
			self.fractionStalled = (self.ribosomeStalls > 0).mean()

		else:
			self.stallingRateTotal = 0
			self.stallingRateMean = 0
			self.stallingRateStd = 0
			self.fractionStalled = 0

	def tableCreate(self, tableWriter):
		pass


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			stallingRateTotal = self.stallingRateTotal,
			stallingRateMean = self.stallingRateMean,
			stallingRateStd = self.stallingRateStd,
			fractionStalled = self.fractionStalled,
			aaCountInSequence = self.aaCountInSequence,
			aaCounts = self.aaCounts,
			actualElongations = self.actualElongations,
			actualElongationHist = self.actualElongationHist,
			elongationsNonTerminatingHist = self.elongationsNonTerminatingHist,
			expectedElongations = self.expectedElongations,
			didTerminate = self.didTerminate,
			didInitialize = self.didInitialize,
			terminationLoss = self.terminationLoss,
			effectiveElongationRate = self.effectiveElongationRate,
			rrn16S_produced = self.rrn16S_produced,
			rrn23S_produced = self.rrn23S_produced,
			rrn5S_produced = self.rrn5S_produced,
			rrn16S_init_prob = self.rrn16S_init_prob,
			rrn23S_init_prob = self.rrn23S_init_prob,
			rrn5S_init_prob = self.rrn5S_init_prob,
			total_rna_init = self.total_rna_init,
			processElongationRate = self.processElongationRate,
			translationSupply = self.translationSupply,
			numTrpATerminated = self.numTrpATerminated,
			probTranslationPerTranscript = self.probTranslationPerTranscript,
			)

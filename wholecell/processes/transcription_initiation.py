#!/usr/bin/env python

"""
Transcription initiation
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class Transcription(wholecell.processes.process.Process):
	""" Transcription """

	# Constructor
	def __init__(self):
		self.meta = {
		"id": "TranscriptionInitiation",
		"name": "TranscriptionInitiation"
		}

		# Constants
		self.rnaPolymeraseTransitionProb = None
		self.promoterBindingProbabilities = None

		super(Transcription, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Transcription, self).initialize(sim, kb)

		## Load constants from Knowledge Base
		self.rnaPolymeraseTransitionProb = kb.rnaPolymeraseTransitionProb
		# kb.rnaPolymeraseTransitionProb = {
		# 	'fromFree' 			: {'toFree' : 0.1, 'toNonspecific': 0.5, 'toSpecific' : 0.4},
		# 	'fromNonspecific'	: {'toFree' : 0.1, 'toNonspecific': 0.5, 'toSpecific' : 0.4},
		# 	'fromSpecific'		: {'toFree' : 0.1, 'toNonspecific': 0.4, 'toSpecific' : 0.4}
		# }


		self.promoterBindingProbabilities = kb.promoterBindingProbabilities
		# kb.promoterBindingProbabilities = numpy.array(len(promoters),
		#	dtype = [('promoterId', 'a'),
		#				'bindingProb', 'f',
		#				'sigmaFactor', 'a'])


		## Create partitions
		# Sigma factors
		self.freeSigma = self.bulkMoleculeView(['RPOD-MONOMER']) # Only Sigma D in here now

		# RNA polymerase
		self.freeRnaPolymerase					= self.bulkMoleculeView(['APORNAP-CPLX[c]'])
		self.nonSpecificallyBoundRnaPolymerase	= self.uniqueMoleculesView('APORNAP-CPLX', {'bindingState' : 'nonSpecific'})
		self.specificallyBoundRnaPolymerase 	= self.uniqueMoleculesView('RNAP70-CPLX'. {'bindingState' : 'specific'})

		# Chromosome
		self.chromosomeAllocation = None
		# Still need to figure out how to create
		# Request RNA polymerase foot print at every promoter
		# Request 


	def calculateRequest(self):
		self.freeSigma.requestAll()
		self.freePolymerase.requestAll()

		# TODO
		self.nonSpecificallyBoundRnaPolymerase
		self.nonSpecificallyBoundRnaPolymerase

		self.chromosomeAllocation


	# Calculate temporal evolution
	def evolveState(self):
		rnaPolymerases = (self.rnapSubunits.counts() // [2, 1, 1, 1]).min()

		ntpEstimate = 1.1 * 4 * self.ntps.counts().min()

		enzLimit = np.min([
			ntpEstimate,
			rnaPolymerases * self.elngRate * self.timeStepSec
			])

		newRnas = 0
		ntpsUsed = np.zeros(4)

		ntpsShape = self.ntps.counts().shape

		rnasCreated = np.zeros_like(self.rnas.counts())

		while enzLimit > 0:
			if not np.any(
					np.all(
						self.ntps.counts() > self.rnaNtCounts,
						axis = 1
						)
					):
				break

			if not np.any(enzLimit > np.sum(self.rnaNtCounts, axis = 1)):
				break

			# If the probabilities of being able to synthesize are sufficiently low, exit the loop
			if np.sum(self.rnaSynthProb[np.all(self.ntps.counts() > self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			if np.sum(self.rnaSynthProb[enzLimit > np.sum(self.rnaNtCounts, axis = 1)]) < 1e-3:
				break

			newIdx = np.where(self.randStream.mnrnd(1, self.rnaSynthProb))[0]

			if np.any(self.ntps.counts() < self.rnaNtCounts[newIdx, :]):
				break

			if enzLimit < np.sum(self.rnaNtCounts[newIdx, :]):
				break

			enzLimit -= np.sum(self.rnaNtCounts[newIdx, :])

			self.ntps.countsDec(
				self.rnaNtCounts[newIdx, :].reshape(ntpsShape)
				)

			self.h2o.countDec(1) # TODO: verify this
			self.ppi.countInc(self.rnaLens[newIdx])
			self.proton.countInc(1) # TODO: verify this

			rnasCreated[newIdx] += 1

		self.rnas.countsInc(rnasCreated)
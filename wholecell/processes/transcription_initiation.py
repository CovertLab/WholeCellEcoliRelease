#!/usr/bin/env python

"""
Transcription initiation
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

from wholecell.utils.package_constants import RNAP_NON_SPECIFICALLY_BOUND_STATE, RNAP_SPECIFICALLY_BOUND_STATE




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
		self.transProb = None

		super(Transcription, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Transcription, self).initialize(sim, kb)

		## Load constants from Knowledge Base
		self.transProb = kb.rnaPolymeraseTransitionProb
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


		self.rnaPolymeraseFootprint = kb.parameters['rnaPolymeraseFootprint']

		## Create partitions
		# Sigma factors
		self.freeSigma = self.bulkMoleculeView(['RPOD-MONOMER']) # Only Sigma D in here now

		# RNA polymerase
		self.freeRnaPolymerase					= self.bulkMoleculeView(['APORNAP-CPLX[c]'])
		self.nonSpecificallyBoundRnaPolymerase	= self.uniqueMoleculesView('APORNAP-CPLX',
			bindingState = ('==', RNAP_NON_SPECIFICALLY_BOUND_STATE))
		self.specificallyBoundRnaPolymerase 	= self.uniqueMoleculesView('RNAP70-CPLX',
			bindingState = ('==', RNAP_SPECIFICALLY_BOUND_STATE))

		# TODO: Active RNAP transitions

		# Chromosome
		self.promoters = self.chromosomeLocationRequest('promoter', self.rnaPolymeraseFootprint)
		self.randomBinding = self.chromosomeRandomRequest(self.rnaPolymeraseFootprint)


	def calculateRequest(self):
		## Compute Markov transitions

		cntSigma = self.freeSigma.total()
		cntFreeRnap = self.freeRnaPolymerase.total()
		cntNonSpecBoundRnap = self.nonSpecificallyBoundRnaPolymerase.total()
		cntSpecBoundRnap = self.specificallyBoundRnaPolymerase.total()
		cntFreePromoters = sum(self.promoters.free())



		expectedTransitionFromFree = cntFreeRnap * (1 - self.transProb['fromFree']['toFree'])
		expectedTransitionFromFreeToSpecific = expectedTransitionFromFree * self.transProb['fromFree']['toSpecific'] / (self.transProb['fromFree']['toSpecific'] + self.transProb['fromFree']['toNonspecific'])


		## Requests
		freeRnapRequest = expectedTransitionFromFree
		sigmaFactorRequest = expectedTransitionFromFreeToSpecific + expectedTransitionFromNonSpecificToSpecific
		if sigmaFactorRequest > cntSigma:
			sigmaFactorRequest = cntSigma

		self.freeRnaPolymerase.requestIs(freeRnapRequest)
		self.freeSigma.requestIs(sigmaFactorRequest)





		self.nonSpecificallyBoundRnaPolymerase.requestIs()
		self.specificallyBoundRnaPolymerase.requestIs()

		self.nonSpecificallyBoundRnaPolymerase
		self.nonSpecificallyBoundRnaPolymerase

		self.chromosomeAllocation


	# Calculate temporal evolution
	def evolveState(self):
		transitions = [self.transProb['fromFree']['toNonspecific'],
						self.transProb['fromFree']['toSpecific'],
						self.transProb['fromFree']['toFree']] # Get rid of free to free prob and normalize to one
		transitions = numpy.cumsum(transitions)

		choices = self.randStream.rand(cntFreeRnap)

		outcomes = len(transitions) - numpy.sum(choices < transitions[:,numpy.newaxis],0) # Index corresponds to transitions

		eachTransition = numpy.bincount(outcomes)

		transitionMatrix = numpy.array([[-1, -1, 0],	# Free
										[1, 0, 0],		# Non-specific
										[0, 1, 0]])		# Specific

		countUpdate = numpy.dot(transitionMatrix, eachTransition)
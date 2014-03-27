#!/usr/bin/env python

"""
Transcription initiation
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

from wholecell.utils.package_constants import RNAP_NON_SPECIFICALLY_BOUND_STATE, RNAP_SPECIFICALLY_BOUND_STATE, RNAP_ACTIVE_STATE

class TranscriptionInitiation(wholecell.processes.process.Process):
	""" TranscriptionInitiation """
	'''
	Transcription initiation is the first step in the cellular processes, which
	produce functional gene products. Transcription initiation begins with the
	recruitment of RNA polymerase (RNAP) to a transcription unit promoter site
	on E. coli's genome with the help of a sigma factor. Once RNAP and sigma
	factor have associated with a promoter site on the DNA active elongation of
	a transcript can begin.
	'''


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
		self.rnaPolymeraseTransitionProb = kb.rnaPolymeraseTransitionProb

		self.promoterBindingProbabilities = kb.promoterBindingProbabilities
		# kb.promoterBindingProbabilities = numpy.array(len(promoters),
		#	dtype = [('promoterId', 'a16'),
		#			('bindingProb','f'),
		#			('sigmaFactor', 'a1'),
		#			('location', 'i'),
		#			('direction','b')])


		self.rnaPolymeraseFootprint = kb.parameters['rnaPolymeraseFootprint']

		## Create partitions
		# Sigma factors
		self.freeSigmaD = self.bulkMoleculeView(['RPOD-MONOMER[c]'])
		self.freeSigmaH = self.bulkMoleculeView(['RPOH-MONOMER[c]'])
		# TODO: Add other sigma factors. Will have to add different sorts of
		# promoters (one for each sigma factor) to the matrix to make this work.

		# RNA polymerase
		self.freeRNAP				= self.bulkMoleculeView(['APORNAP-CPLX[c]'])
		self.nonSpecificRNAP		= self.uniqueMoleculesView('APORNAP-CPLX', bindingState = ('==', RNAP_NON_SPECIFICALLY_BOUND_STATE))
		self.specificRNAP_sigmaD 	= self.uniqueMoleculesView('RNAP70-CPLX', bindingState = ('==', RNAP_SPECIFICALLY_BOUND_STATE))
		self.specificRNAP_sigmaH 	= self.uniqueMoleculesView('RNAP32-CPLX',bindingState = ('==', RNAP_SPECIFICALLY_BOUND_STATE))
		self.activeRNAP_sigmaD		= self.uniqueMoleculesView('RNAP70-CPLX', bindingState = ('==', RNAP_ACTIVE_BOUND_STATE))
		self.activeRNAP_sigmaH 		= self.uniqueMoleculesView('RNAP32-CPLX', bindingState = ('==', RNAP_ACTIVE_BOUND_STATE))

		# Chromosome
		self.promoters_sigmaD = self.chromosomeLocationRequest('promoter', self.rnaPolymeraseFootprint)
		self.promoters_sigmaH
		self.randomBinding = self.chromosomeRandomRequest(self.rnaPolymeraseFootprint)

		# Transcripts
		# TODO: Create container views for transcripts? Or at least a pointer
		# to the transcript state so that they can be created?


	def calculateRequest(self):
		'''
		calculateRequest

		Requests are calculated based on expected transitions in the Markov model
		representing the different states that an RNA polymerase can exist in: 
		free, non-specifically bound, specifically bound, or actively transcribing.
		In this way the process only requests what will be used in the next
		time step.
		'''

		## Compute Markov transitions to calculate requests
		###################################################
		# n = number of states
		# v = flux through Markov transitions (dimension (n-1,1))
		# P = probability matrix for Markov transitions with columns
		#	  summing to unity (dimension (n,n))
		# Q = probability matrix with each reaction having its own row.
		#	  Just a transformation of P. (dimension (n^n,n))
		# m = counts of Markov states (dimension (n,1))
		# S = stoichiometric matrix for Markov transitions (dimension (n, n^n))
		# 
		# Calculate fluxes through transitions
		# 	v = Qm
		# Calculate change in states
		# 	m(t+1) - m(t) = Sv
		# Substitute
		# 	m(t+1) - m(t) = SQm(t)
		# Requests are things that are needed by the process not produced
		# so we need the negative max relative to zero.
		# 	req = max(0, -(m(t+1) - m(t))) = max(0, -SQm(t))
		# Which in Python is written
		#	numpy.fmax(0, -1. * numpy.dot(numpy.dot(S,Q),m))

		# Index/Species:
		# 	0	/	free RNAP
		#	1 	/	non specifically bound RNAP
		#	2 	/	specifically bound RNAP
		#	3	/	activly transcribing
		# --------------------- below only in S matrix
		#	4	/	sigma factor
		#	5	/	free promoters

		# TODO: Load from KB
		F_idx = 0
		NS_idx = 1
		S_idx = 2
		A_idx = 3
		sigma_idx = 4
		promoter_idx = 5

		m = numpy.array([self.freeRNAP.total(),
						self.nonSpecificRNAP.total(),
						self.specificRNAP_sigmaD.total() + self.specificRNAP_sigmaH,
						self.activeRNAP_sigmaH.total() + self.activeRNAP_sigmaH
						])
		
		# TODO: Load from KB
		# TODO: Shouldn't have a transition from A --> F that should
		# go into terminatoin.
		#			From 	F 		NS 		S 	A
		P = numpy.array([[	0.1,	0.1,	0,	0.1	], # to F
						 [	0.9,	0.4,	0,	0	], # to NS
						 [	0,		0.5,	0,	0	], # to S
						 [	0,		0,		1,	0.9	]])# to A

		P_shape = P.shape[0]
		Q = numpy.zeros((P_shape**2,P_shape))
		for i in range(P_shape):
			Q[i*P_shape:(i+1)*P_shape, i] = P[:,i]

		# TODO: Load from KB
		# Assuming that specifically and non-specifically bound RNAP when staying in same state
		# don't move.
		#					F 	F 	 F   F  NS  NS  NS  NS   S   S   S   S   A   A   A   A
		#					to
		#					F 	NS 	 S   A   F  NS   S   A   F  NS   S   A   F  NS   S   A
		S = numpy.array([   [0, -1, -1, -1,  1,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0],	# Free
							[0,  1,  0,  0, -1,  0, -1, -1,  0,  1,  0,  0,  0,  1,  0,  0],	# Non-specifically bound
							[0,  0,  1,  0,  0,  0,  1,  0, -1, -1,  0, -1,  0,  0,  1,  0],	# Specifically bound
							[0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  1, -1, -1, -1,  0],	# Active
							[0,  0, -1, -1,  0,  0, -1, -1,  1,  1,  0,  0,  1,  0,  0,  0],	# Sigma
							[0,  0, -1,  0,  0,  0, -1,  0,  1,  1,  0,  1,  0,  1, -1,  0]])	# Promoter

		speciesUsedInTransitions = self.randStream.stochasticRound((numpy.dot(numpy.dot(S,Q),m)))


		## Limit requests by available sigma factors and promoters
		##########################################################
		sigmaAndPromoterIdx = numpy.array([sigma_idx, promoter_idx])
		if numpy.min(speciesUsedInTransitions[sigmaAndPromoterIdx]) < -1 * numpy.min(self.freeSigmaD.total(), len(self.promoters.free())):
			# Limit NS --> S transition reaction
			# NS + Sigma + Promoter --> S
			desiredAmount = numpy.min(speciesUsedInTransitions[sigmaAndPromoterIdx])
			maxAvailable = -1 * numpy.min([self.freeSigmaD.total(), len(self.promoters.free()])
			difference = maxAvailable - desiredAmount # Will be positive
			speciesUsedInTransitions[S_idx] 		-= difference
			speciesUsedInTransitions[NS_idx] 		+= difference
			speciesUsedInTransitions[sigma_idx] 	+= difference
			speciesUsedInTransitions[promoter_idx]	+= difference

		requests = numpy.fmax(0, -1. * speciesUsedInTransitions)

		## Set requests
		###############

		self.freeRNAP.requestIs(requests[F_idx])
		self.nonSpecificRNAP.requestIs(requests[NS_idx])

		# TODO: Scale here somehow? So that the total request is proprotional
		# to the amount of each type of sigma that is required for promoters?
		self.specificRNAP_sigmaD.requestIs(requests[S_idx])
		self.specificRNAP_sigmaH.requestIs()


		# TODO: Should probably just remove this because this process will
		# only really be creating A not using it for anything. For now, it is
		# here.
		self.activeRNAP.requestIs(requests[A_idx])

		# TODO: Write chromsome allocation for promoters and random binding

	# Calculate temporal evolution
	def evolveState(self):
		'''
		evolveState

		Expected transitions in the Markov model have already been calculated
		in the calculateRequest function...
		'''

		## Initiate specifically bound polymerases
		##########################################

		for molecule in self.specificRNAP:
			molecule.attrIs(bindingState = RNAP_ACTIVE_STATE)
			# TODO: Create transcript
			# TODO: Set RNAP attributes so that it is associated with that transcript.


		## Specifically bound polymerases becoming non-specitically bound/free
		######################################################################
		# NOTE: Not in the model right now but could/should be


		## Free or non-specifically bound polymerases becoming specifically bound
		#########################################################################






# NOTES
# RNA polymerases - unique attributes:
# - bindingState					- (F, NS, S, A, T)
# - associatedTranscript 			- Idx of transcript it is elongating
# - associatedTranscriptionUnit		- Idx of transcription unit it is elongating
# - associatedSigmaFactor			- Sigma factor it is associated with
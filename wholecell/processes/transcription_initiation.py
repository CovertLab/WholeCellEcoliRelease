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
		self.rnaPolymeraseTransitionProb = kb.rnaPolymeraseTransitionProb

		self.promoterBindingProbabilities = kb.promoterBindingProbabilities
		# kb.promoterBindingProbabilities = numpy.array(len(promoters),
		#	dtype = [('promoterId', 'a'),
		#				'bindingProb', 'f',
		#				'sigmaFactor', 'a'])


		self.rnaPolymeraseFootprint = kb.parameters['rnaPolymeraseFootprint']

		## Create partitions
		# Sigma factors
		self.freeSigma = self.bulkMoleculeView(['RPOD-MONOMER']) # Only Sigma D in here now
		# TODO: Add other initiation factors, or factors required to go from,
		# specifically bound to active.

		# RNA polymerase
		self.freeRnaPolymerase					= self.bulkMoleculeView(['APORNAP-CPLX[c]'])
		self.nonSpecificallyBoundRnaPolymerase	= self.uniqueMoleculesView('APORNAP-CPLX',
			bindingState = ('==', RNAP_NON_SPECIFICALLY_BOUND_STATE))
		self.specificallyBoundRnaPolymerase 	= self.uniqueMoleculesView('RNAP70-CPLX',
			bindingState = ('==', RNAP_SPECIFICALLY_BOUND_STATE))

		self.activeRnaPolymerase 				= self.uniqueMoleculesView('RNAP70-CPLX',
			bindingState = ('==', RNAP_ACTIVE_BOUND_STATE))

		# Chromosome
		self.promoters = self.chromosomeLocationRequest('promoter', self.rnaPolymeraseFootprint)
		self.randomBinding = self.chromosomeRandomRequest(self.rnaPolymeraseFootprint)


	def calculateRequest(self):
		## Compute Markov transitions to calculate requests
		# 
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

		m = numpy.array([self.freeRnaPolymerase.total(),
						self.nonSpecificallyBoundRnaPolymerase.total(),
						self.specificallyBoundRnaPolymerase.total(),
						self.activeRnaPolymerase.total()
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

		sigmaAndPromoterIdx = numpy.array([sigma_idx, promoter_idx])
		if numpy.min(speciesUsedInTransitions[sigmaAndPromoterIdx]) < -1 * numpy.min(self.freeSigma.total(), len(self.promoters.free())):
			# Limit NS --> S transition reaction
			# NS + Sigma + Promoter --> S
			desiredAmount = numpy.min(speciesUsedInTransitions[sigmaAndPromoterIdx])
			maxAvailable = -1 * numpy.min([self.freeSigma.total(), len(self.promoters.free()])
			difference = maxAvailable - desiredAmount # Will be positive
			speciesUsedInTransitions[S_idx] 		-= difference
			speciesUsedInTransitions[NS_idx] 		+= difference
			speciesUsedInTransitions[sigma_idx] 	+= difference
			speciesUsedInTransitions[promoter_idx]	+= difference

		requests = numpy.fmax(0, -1. * speciesUsedInTransitions)

		self.freeRnaPolymerase.requestIs(requests[F_idx])
		self.nonSpecificallyBoundRnaPolymerase.requestIs(requests[NS_idx])
		self.specificallyBoundRnaPolymerase.requestIs(requests[S_idx])
		# TODO: Should probably just remove this because this process will
		# only really be creating A not using it for anything. For now, it is
		# here.
		self.activeRnaPolymerase.requestIs(requests[A_idx])

		self.chromosomeAllocation

	# Calculate temporal evolution
	def evolveState(self):
		pass
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
	def initialize(self, sim, kb):
		super(TranscriptInitiation, self).initialize(sim, kb)

		# Load parameters

		self.rnaSynthProb = kb.rnaData['synthProb']

		self.activationProb = kb.transcriptionActivationRate.asNumber(1/units.s) * self.timeStepSec # TODO: consider the validity of this math

		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')

		self.inactiveRnaPolys = self.bulkMoleculeView("APORNAP-CPLX[c]")


	def calculateRequest(self):
		self.inactiveRnaPolys.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# Sample a multinomial distribution of synthesis probabilities to 
		# determine what molecules are initialized

		inactiveRnaPolyCount = self.inactiveRnaPolys.count()

		rnaPolyToActivate = np.int64(self.activationProb * inactiveRnaPolyCount)

		if rnaPolyToActivate == 0:
			return

		nNewRnas = self.randomState.multinomial(rnaPolyToActivate,
			self.rnaSynthProb)

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

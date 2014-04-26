#!/usr/bin/env python

"""
UniqueTranscriptInitiation

Transcription initiation sub-model.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/26/14
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

import itertools

class UniqueTranscriptInitiation(wholecell.processes.process.Process):
	""" UniqueTranscriptInitiation """

	_name = "UniqueTranscriptInitiation"

	# Constructor
	def __init__(self):
		# Constants
		self.rnaIds = None
		self.rnaNtCounts = None
		self.rnaSynthProb = None

		super(UniqueTranscriptInitiation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(UniqueTranscriptInitiation, self).initialize(sim, kb)

		# Load parameters

		# TODO: evaluate the roles of enzymes in initiation

		# enzIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]

		self.rnaIds = kb.rnaData['id']
		self.rnaNtCounts = kb.rnaData['countsAUCG']
		self.rnaSynthProb = kb.rnaData['synthProb']

		self.initiationRate = 100 # TODO: move parameter to KB and fit!

		# Views

		self.rnas = self.uniqueMoleculesView('rnaTranscript')

		# self.rnapSubunits = self.bulkMoleculesView(enzIds)


	def calculateRequest(self):
		# No request, since we're only creating 'empty' transcripts
		pass


	# Calculate temporal evolution
	def evolveState(self):
		# Sample a multinomial distribution of synthesis probabilities to 
		# determine what molecules are initializaed

		nNewRnas = self.randStream.mnrnd(self.initiationRate,
			self.rnaSynthProb)

		nonzeroCount = (nNewRnas > 0)

		for nNew, rnaId, ntCounts in itertools.izip(
				nNewRnas[nonzeroCount],
				self.rnaIds[nonzeroCount],
				self.rnaNtCounts[nonzeroCount]
				):

			self.rnas.moleculesNew(
				'rnaTranscript', nNew,
				rnaId = rnaId,
				requiredAUCG = ntCounts
				)

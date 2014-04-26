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

		enzIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]

		self.rnaNtCounts = kb.rnaData['countsAUCG']
		self.rnaSynthProb = kb.rnaData['synthProb']

		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')

		self.rnapSubunits = self.bulkMoleculesView(enzIds)


	def calculateRequest(self):
		self.rnapSubunits.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# Sample a multinomial distribution of synthesis probabilities to 
		# determine what molecules are initialized

		inactiveRnaPolys = (self.rnapSubunits.counts() // [2, 1, 1, 1]).min()

		nNewRnas = self.randStream.mnrnd(inactiveRnaPolys,
			self.rnaSynthProb)

		# Create the active RNA polymerases

		nonzeroCount = (nNewRnas > 0)

		for rnaIndex, (nNew, ntCounts) in enumerate(itertools.izip(
				nNewRnas[nonzeroCount],
				self.rnaNtCounts[nonzeroCount]
				)):

			self.activeRnaPolys.moleculesNew(
				'activeRnaPoly', nNew,
				rnaIndex = rnaIndex,
				requiredAUCG = ntCounts
				)

		self.rnapSubunits.countsIs(0)

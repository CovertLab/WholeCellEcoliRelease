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
		# Parameters
		self.rnaSynthProb = None

		# Views
		self.activeRnaPolys = None
		self.rnapSubunits = None

		super(UniqueTranscriptInitiation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(UniqueTranscriptInitiation, self).initialize(sim, kb)

		# Load parameters

		enzIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]

		self.rnaSynthProb = kb.rnaData['synthProb'].to('dimensionless').magnitude

		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')

		self.rnapSubunits = self.bulkMoleculesView(enzIds)


	def calculateRequest(self):
		self.rnapSubunits.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# Sample a multinomial distribution of synthesis probabilities to 
		# determine what molecules are initialized

		inactiveRnaPolyCount = (self.rnapSubunits.counts() // [2, 1, 1, 1]).min()

		if inactiveRnaPolyCount == 0:
			return

		nNewRnas = self.randomState.multinomial(inactiveRnaPolyCount,
			self.rnaSynthProb)

		nonzeroCount = (nNewRnas > 0)

		assert nNewRnas.sum() == inactiveRnaPolyCount

		# Build list of RNA indexes

		rnaIndexes = np.empty(inactiveRnaPolyCount, np.int64)

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
			inactiveRnaPolyCount
			)

		activeRnaPolys.attrIs(
			rnaIndex = rnaIndexes
			)

		self.rnapSubunits.countsDec(
			nNewRnas.sum() * np.array([2, 1, 1, 1], np.int)
			)

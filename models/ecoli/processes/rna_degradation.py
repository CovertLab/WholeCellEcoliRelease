#!/usr/bin/env python

"""
RnaDegradation

RNA degradation sub-model. Encodes molecular simulation of RNA degradation as a Poisson process

TODO:
- handle complexes

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.constants import REQUEST_PRIORITY_DEGRADATION
from wholecell.utils import units

class RnaDegradation(wholecell.processes.process.Process):
	""" RnaDegradation """

	_name = "RnaDegradation"

	# Constructor
	def __init__(self):
		# Parameters
		self.rnaLens = None			# RNA lengths
		self.rnaDegRates = None		# RNA degradation rates (1/s)
		self.rnaDegSMat = None		# RNA degradation stoichiometry matrix [metabolite x rna]

		# Views
		self.metabolites = None
		self.rnas = None
		self.rnase = None

		super(RnaDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(RnaDegradation, self).initialize(sim, kb)

		metaboliteIds = ["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]",
			"H2O[c]", "PPI[c]", "H[c]"]

		nmpIdxs = np.arange(0, 4)
		h2oIdx = metaboliteIds.index('H2O[c]')
		ppiIdx = metaboliteIds.index('PPI[c]')
		hIdx = metaboliteIds.index('H[c]')

		rnaIds = kb.process.transcription.rnaData['id']

		# Rna
		self.rnaDegRates = kb.process.transcription.rnaData['degRate'].asNumber()
		self.rnaLens = kb.process.transcription.rnaData['length'].asNumber()

		# Build stoichiometric matrix
		# TODO: account for NTP on 5' end
		self.rnaDegSMat = np.zeros((len(metaboliteIds), len(rnaIds)), np.int64)
		self.rnaDegSMat[nmpIdxs, :] = units.transpose(kb.process.transcription.rnaData['countsACGU']).asNumber()
		# self.rnaDegSMat[h2oIdx, :]  = -(self.rnaLens - 1)
		self.rnaDegSMat[h2oIdx, :]  = -self.rnaLens # using one additional water to hydrolyze PPI on 5' end
		self.rnaDegSMat[ppiIdx, :]    =  1
		self.rnaDegSMat[hIdx, :] = self.rnaLens
		
		# Views
		self.metabolites = self.bulkMoleculesView(metaboliteIds)		
		self.rnas = self.bulkMoleculesView(rnaIds)
		self.rnase = self.bulkMoleculeView('EG11259-MONOMER[c]')

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)


	# Calculate temporal evolution

	def calculateRequest(self):
		nRNAsToDegrade = np.fmin(
			self.randomState.poisson(self.rnaDegRates * self.rnas.total() * self.timeStepSec),
			self.rnas.total()
			)

		# nReactions = np.dot(self.rnaLens, nRNAsToDegrade)

		# self.h2o.requestIs(nReactions)
		self.rnas.requestIs(nRNAsToDegrade)
		self.rnase.requestAll()

		metaboliteUsage = np.fmax(
			-np.dot(self.rnaDegSMat, nRNAsToDegrade),
			0
			)

		self.metabolites.requestIs(metaboliteUsage)
		

	def evolveState(self):
		# Check if RNAse R expressed
		if self.rnase.count() == 0:
			return

		# Degrade RNA
		self.metabolites.countsInc(np.dot(
			self.rnaDegSMat,
			self.rnas.counts()
			))

		self.rnas.countsIs(0)

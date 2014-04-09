#!/usr/bin/env python

"""
RnaDegradation

RNA degradation sub-model. Encodes molecular simulation of RNA degradation as a Poisson process

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

import numpy as np

import wholecell.processes.process


class RnaDegradation(wholecell.processes.process.Process):
	""" RnaDegradation """

	_metaboliteIds = None
	_rnaIds = None

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "RnaDegradation",
			"name": "RNA degradation"
		}

		# Constants
		self.rnaLens = None			# RNA lengths
		self.rnaDegRates = None		# RNA degradation rates (1/s)
		self.rnaDegSMat = None		# RNA degradation stoichiometry matrix [metabolite x rna]

		super(RnaDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb, kb2):
		super(RnaDegradation, self).initialize(sim, kb, kb2)

		self._metaboliteIds = ["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]",
			"H2O[c]", "H[c]", "ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"]

		self._nmpIdxs = np.arange(0, 4)
		self._h2oIdx = self._metaboliteIds.index('H2O[c]')
		self._hIdx = self._metaboliteIds.index('H[c]')

		self._rnaIds = [x["id"] + "[c]" for x in kb.rnas]

		# Rna
		self.rnaDegRates = np.log(2) / np.array([x["halfLife"] for x in kb.rnas])

		self.rnaLens = np.sum(np.array([x["ntCount"] for x in kb.rnas]), axis = 1)

		self.rnaDegSMat = np.zeros((len(self._metaboliteIds), len(self._rnaIds)), np.int64)
		self.rnaDegSMat[self._nmpIdxs, :] = np.transpose(np.array([x["ntCount"] for x in kb.rnas]))
		self.rnaDegSMat[self._h2oIdx, :]  = -(np.sum(self.rnaDegSMat[self._nmpIdxs, :], axis = 0) - 1)
		self.rnaDegSMat[self._hIdx, :]    =  (np.sum(self.rnaDegSMat[self._nmpIdxs, :], axis = 0) - 1)

		# Views
		self.metabolites = self.bulkMoleculesView(self._metaboliteIds)

		self.nmps = self.bulkMoleculesView(["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"])
		self.h2o = self.bulkMoleculeView('H2O[c]')
		self.proton = self.bulkMoleculeView('H[c]')

		self.rnas = self.bulkMoleculesView(self._rnaIds)

		self.rnase = self.bulkMoleculeView('EG11259-MONOMER[c]')


	# Calculate temporal evolution

	def calculateRequest(self):
		nRNAsToDegrade = np.fmin(
			self.randStream.poissrnd(self.rnaDegRates * self.rnas.total() * self.timeStepSec),
			self.rnas.total()
			)

		nReactions = np.dot(self.rnaLens, nRNAsToDegrade)

		self.h2o.requestIs(nReactions)
		self.rnas.requestIs(nRNAsToDegrade)
		self.rnase.requestAll()
		

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

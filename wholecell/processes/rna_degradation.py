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
	def initialize(self, sim, kb):
		super(RnaDegradation, self).initialize(sim, kb)

		self._metaboliteIds = ["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]",
			"H2O[c]", "H[c]", "ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"]

		# self._ntpIdxs = np.arange(6, 10)
		self._nmpIdxs = np.arange(0, 4)
		self._h2oIdx = self._metaboliteIds.index('H2O[c]')
		self._hIdx = self._metaboliteIds.index('H[c]')

		self._rnaIds = [x["id"] + "[c]" for x in kb.rnas]

		mc = sim.states['BulkMolecules']

		self.bulkMoleculesPartition.initialize(self._metaboliteIds + self._rnaIds + ["EG11259-MONOMER[c]"])

		# Metabolites
		self.bulkMoleculesPartition.metabolites = self.bulkMoleculesPartition.countsView(self._metaboliteIds)

		self.bulkMoleculesPartition.nmps = self.bulkMoleculesPartition.countsView(["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"])
		self.bulkMoleculesPartition.h2oMol = self.bulkMoleculesPartition.countView('H2O[c]')
		self.bulkMoleculesPartition.hMol = self.bulkMoleculesPartition.countView('H[c]')

		# Rna
		self.bulkMoleculesPartition.rnas = self.bulkMoleculesPartition.countsView(self._rnaIds)

		self.rnaView = mc.countsView(self._rnaIds)

		self.rnaDegRates = np.log(2) / np.array([x["halfLife"] for x in kb.rnas])

		self.rnaLens = np.sum(np.array([x["ntCount"] for x in kb.rnas]), axis = 1)

		self.rnaDegSMat = np.zeros((len(self._metaboliteIds), len(self._rnaIds)))
		self.rnaDegSMat[self._nmpIdxs, :] = np.transpose(np.array([x["ntCount"] for x in kb.rnas]))
		self.rnaDegSMat[self._h2oIdx, :]  = -(np.sum(self.rnaDegSMat[self._nmpIdxs, :], axis = 0) - 1)
		self.rnaDegSMat[self._hIdx, :]    =  (np.sum(self.rnaDegSMat[self._nmpIdxs, :], axis = 0) - 1)

		# Proteins
		self.bulkMoleculesPartition.rnaseRMol = self.bulkMoleculesPartition.countView('EG11259-MONOMER[c]')

		# Views
		self.metabolites = self.bulkMoleculesView(self._metaboliteIds)
		self.nmps = self.bulkMoleculesView(["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"])
		self.h2o = self.bulkMoleculeView('H2O[c]')
		self.proton = self.bulkMoleculeView('H[c]')
		self.rnas = self.bulkMoleculesView(self._rnaIds)
		self.rnase = self.bulkMoleculeView('EG11259-MONOMER[c]')


	def requestBulkMolecules(self):
		self.bulkMoleculesPartition.h2oMol.countIs(
			np.dot(self.rnaLens, self.rnaDegRates * self.rnaView.counts()) * self.timeStepSec
			)
		
		self.bulkMoleculesPartition.rnas.countsIs(
			self.randStream.poissrnd(self.rnaDegRates * self.rnaView.counts() * self.timeStepSec)
			)

		self.bulkMoleculesPartition.rnaseRMol.countInc(1)


	def calculateRequest(self):
		nRNAsToDegrade = np.fmin(
			self.randStream.poissrnd(self.rnaDegRates * self.rnas.total() * self.timeStepSec),
			self.rnas.total()
			)

		nReactions = np.dot(self.rnaLens, nRNAsToDegrade)

		self.h2o.requestIs(nReactions)
		self.rnas.requestIs(nRNAsToDegrade)
		self.rnase.requestAll()
		

	# Calculate temporal evolution
	def evolveState(self):
		# Check if RNAse R expressed
		if self.bulkMoleculesPartition.rnaseRMol.count() == 0:
			return

		# Degrade RNA
		self.bulkMoleculesPartition.metabolites.countsInc(
			np.dot(self.rnaDegSMat, self.bulkMoleculesPartition.rnas.counts())
			)

		self.bulkMoleculesPartition.rnas.countsIs(0)

		# print "NTP recycling: %s" % str(self.metabolite.counts[self.metabolite.idx["ntps"]])

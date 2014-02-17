#!/usr/bin/env python

"""
RnaDegradation

RNA degradation sub-model. Encodes molecular simulation of RNA degradation as a Poisson process

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

import numpy

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

		# self._ntpIdxs = numpy.arange(6, 10)
		self._nmpIdxs = numpy.arange(0, 4)
		self._h2oIdx = self._metaboliteIds.index('H2O[c]')
		self._hIdx = self._metaboliteIds.index('H[c]')

		self._rnaIds = [x["id"] + ":nascent[c]" for x in kb.rnas] + [x["id"] + "[c]" for x in kb.rnas]

		mc = sim.states['MoleculeCounts']

		mc.setPartition(self, self._metaboliteIds + self._rnaIds + ["EG11259-MONOMER[c]"])

		# Metabolites
		self.mcPartition.metabolites = self.mcPartition.countsBulkViewNew(self._metaboliteIds)

		self.mcPartition.nmps = self.mcPartition.countsBulkViewNew(["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"])
		self.mcPartition.h2oMol = self.mcPartition.molecule('H2O[c]')
		self.mcPartition.hMol = self.mcPartition.molecule('H[c]')

		# Rna
		self.mcPartition.rnas = self.mcPartition.countsBulkViewNew(self._rnaIds)

		self.rnaView = mc.countsBulkViewNew(self._rnaIds)

		self.rnaDegRates = numpy.log(2) / numpy.array([x["halfLife"] for x in kb.rnas] * 2)

		self.rnaLens = numpy.sum(numpy.array([x["ntCount"] for x in kb.rnas] * 2), axis = 1)

		self.rnaDegSMat = numpy.zeros((len(self._metaboliteIds), len(self._rnaIds)))
		self.rnaDegSMat[self._nmpIdxs, :] = numpy.transpose(numpy.array([x["ntCount"] for x in kb.rnas] * 2))
		self.rnaDegSMat[self._h2oIdx, :]  = -(numpy.sum(self.rnaDegSMat[self._nmpIdxs, :], axis = 0) - 1)
		self.rnaDegSMat[self._hIdx, :]    =  (numpy.sum(self.rnaDegSMat[self._nmpIdxs, :], axis = 0) - 1)

		# Proteins
		self.mcPartition.rnaseRMol = self.mcPartition.molecule('EG11259-MONOMER[c]')


	# Calculate temporal evolution
	def evolveState(self):
		# Check if RNAse R expressed
		if self.mcPartition.rnaseRMol.countBulk() == 0:
			return

		# Degrade RNA
		self.mcPartition.metabolites.countsBulkInc(
			numpy.dot(self.rnaDegSMat, self.mcPartition.rnas.countsBulk())
			)

		self.mcPartition.rnas.countsBulkIs(0)

		# print "NTP recycling: %s" % str(self.metabolite.counts[self.metabolite.idx["ntps"]])


	def requestMoleculeCounts(self):
		self.mcPartition.h2oMol.countBulkIs(
			numpy.dot(self.rnaLens, self.rnaDegRates * self.rnaView.countsBulk()) * self.timeStepSec
			)
		
		self.mcPartition.rnas.countsBulkIs(
			self.randStream.poissrnd(self.rnaDegRates * self.rnaView.countsBulk() * self.timeStepSec)
			)

		self.mcPartition.rnaseRMol.countBulkInc(1)

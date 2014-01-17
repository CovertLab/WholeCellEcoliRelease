#!/usr/bin/env python

"""
RnaDegradation

RNA degradation sub-model. Encodes molecular simulation of RNA degradation as a Poisson process

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

import numpy

import wholecell.sim.process.Process

class RnaDegradation(wholecell.sim.process.Process.Process):
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

		self._ntpIdxs = numpy.arange(6, 10)
		self._nmpIdxs = numpy.arange(0, 4)
		self._h2oIdx = self._metaboliteIds.index('H2O[c]')
		self._hIdx = self._metaboliteIds.index('H[c]')

		self._rnaIds = [x["id"] + ":nascent[c]" for x in kb.rnas] + [x["id"] + ":mature[c]" for x in kb.rnas]

		mc = sim.getState('MoleculeCounts')

		# Metabolites
		self.metabolite = mc.addPartition(self, self._metaboliteIds, self.calcReqMetabolites)

		self.nmpView = self.metabolite.countsBulkViewNew(["AMP", "CMP", "GMP", "UMP"])
		self.ntpView = self.metabolite.countsBulkViewNew(["ATP", "CTP", "GTP", "UTP"])
		self.h2oMol = self.metabolite.molecule('H2O:mature', 'merged') # TODO: fix compartment referencing in partitions
		self.hMol = self.metabolite.molecule('H:mature', 'merged') # TODO: fix compartment referencing in partitions

		# Rna
		self.rna = mc.addPartition(self, self._rnaIds ,self.calcReqRna, True)

		self.rnaView = mc.countsBulkViewNew(self._rnaIds)

		self.rnaDegRates = numpy.log(2) / numpy.array([x["halfLife"] for x in kb.rnas] * 2)

		self.rnaLens = numpy.sum(numpy.array([x["ntCount"] for x in kb.rnas] * 2), axis = 1)

		self.rnaDegSMat = numpy.zeros((len(self._metaboliteIds), len(self._rnaIds)))
		self.rnaDegSMat[self._ntpIdxs, :] = numpy.transpose(numpy.array([x["ntCount"] for x in kb.rnas] * 2))
		self.rnaDegSMat[self._h2oIdx, :]  = -(numpy.sum(self.rnaDegSMat[self._nmpIdxs, :], axis = 0) - 1)
		self.rnaDegSMat[self._hIdx, :]    =  (numpy.sum(self.rnaDegSMat[self._nmpIdxs, :], axis = 0) - 1)

		# Proteins
		self.enzyme = mc.addPartition(self, ["EG11259-MONOMER:mature[c]"], self.calcReqEnzyme)

		self.rnaseRMol = self.enzyme.molecule('EG11259-MONOMER:mature', 'merged')


	# Calculate temporal evolution
	def evolveState(self):
		# Check if RNAse R expressed
		if self.rnaseRMol.countBulk() == 0:
			return

		# Degrade RNA
		self.metabolite.countsBulkInc(
			numpy.dot(self.rnaDegSMat, self.rnaView.countsBulk())
			)

		self.rna.countsBulkIs(0)

		# print "NTP recycling: %s" % str(self.metabolite.counts[self.metabolite.idx["ntps"]])

		
	# Calculate needed metabolites
	def calcReqMetabolites(self, request):
		request.countsBulkIs(0)

		self.h2oMol.countBulkIs(
			numpy.dot(self.rnaLens, self.rnaDegRates * self.rnaView.countsBulk()) * self.timeStepSec
			)


	# Calculate needed RNA
	def calcReqRna(self, request):
		request.countsBulkIs(
			self.randStream.poissrnd(self.rnaDegRates * self.rnaView.countsBulk() * self.timeStepSec)
			)


	# Calculate needed proteins
	def calcReqEnzyme(self, request):
		request.countsBulkIs(1)

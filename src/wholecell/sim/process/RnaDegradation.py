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

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "RnaDegradation",
			"name": "RNA degradation"
		}

		# References to states
		self.metabolite = None
		self.rna = None
		self.enzyme = None

		# Constants
		self.rnaLens = None			# RNA lengths
		self.rnaDegRates = None		# RNA degradation rates (1/s)
		self.rnaDegSMat = None		# RNA degradation stoichiometry matrix [metabolite x rna]

		super(RnaDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(RnaDegradation, self).initialize(sim, kb)

		# Metabolites
		self.metabolite = sim.getState("MoleculeCounts").addPartition(self, [
			"AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]", "H2O[c]", "H[c]",
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"
			], self.calcReqMetabolites)
		self.metabolite.idx["nmps"] = self.metabolite.getIndex(["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"])[0]
		self.metabolite.idx["ntps"] = self.metabolite.getIndex(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])[0]
		self.metabolite.idx["h2o"] = self.metabolite.getIndex("H2O[c]")[0]
		self.metabolite.idx["h"] = self.metabolite.getIndex("H[c]")[0]

		# Rna
		self.rna = sim.getState("MoleculeCounts").addPartition(self, 
			[x["id"] + ":nascent[c]" for x in kb.rnas] + [x["id"] + ":mature[c]" for x in kb.rnas]
			, self.calcReqRna, True)

		self.rnaDegRates = numpy.log(2) / numpy.array([x["halfLife"] for x in kb.rnas] * 2)

		self.rnaLens = numpy.sum(numpy.array([x["ntCount"] for x in kb.rnas] * 2), axis = 1)

		self.rnaDegSMat = numpy.zeros((len(self.metabolite.ids), len(self.rna.ids)))
		self.rnaDegSMat[self.metabolite.idx["ntps"], :] = numpy.transpose(numpy.array([x["ntCount"] for x in kb.rnas] * 2))
		self.rnaDegSMat[self.metabolite.idx["h2o"], :]  = -(numpy.sum(self.rnaDegSMat[self.metabolite.idx["nmps"], :], axis = 0) - 1)
		self.rnaDegSMat[self.metabolite.idx["h"], :]    =  (numpy.sum(self.rnaDegSMat[self.metabolite.idx["nmps"], :], axis = 0) - 1)

		# Proteins
		self.enzyme = sim.getState("MoleculeCounts").addPartition(self, ["EG11259-MONOMER:mature[c]"], self.calcReqEnzyme)
		self.enzyme.idx["rnaseR"] = self.enzyme.getIndex(["EG11259-MONOMER:mature[c]"])[0]

	# Calculate needed metabolites
	def calcReqMetabolites(self):
		val = numpy.zeros(self.metabolite.fullCounts.shape)
		val[self.metabolite.idx["h2o"]] = numpy.dot(self.rnaLens, self.rnaDegRates * self.rna.fullCounts) * self.timeStepSec
		return val

	# Calculate needed RNA
	def calcReqRna(self):
		return self.randStream.poissrnd(self.rnaDegRates * self.rna.fullCounts * self.timeStepSec)

	# Calculate needed proteins
	def calcReqEnzyme(self):
		return numpy.ones(self.enzyme.fullCounts.shape)

	# Calculate temporal evolution
	def evolveState(self):
		# Check if RNAse R expressed
		if self.enzyme.counts[self.enzyme.idx["rnaseR"]] == 0:
			return

		# Degrade RNA
		self.metabolite.counts += numpy.dot(self.rnaDegSMat, self.rna.counts)
		self.rna.counts[:] = 0

		print "NTP recycling: %s" % str(self.metabolite.counts[self.metabolite.idx["ntps"]])
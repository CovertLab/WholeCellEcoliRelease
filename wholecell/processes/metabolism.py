#!/usr/bin/env python

"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

from wholecell.reconstruction.fitter import normalize, countsFromMassAndExpression
from wholecell.utils.random import stochasticRound

class Metabolism(wholecell.processes.process.Process):
	""" Metabolism """

	_name = "Metabolism"

	# Constructor
	def __init__(self):

		# Parameters
		self.nAvogadro = None
		self.initialDryMass = None
		self.cellCycleLen = None
		self.wildtypeBiomassReaction = None
		self.wildtypeIds = None

		# Views
		self.biomassMetabolites = None
		self.ppi = None
		self.ntpsdntps = None
		self.ntps = None
		self.nmps = None

		super(Metabolism, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Metabolism, self).initialize(sim, kb)

		# Load constants
		self.nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
		self.initialDryMass = kb.avgCellDryMassInit.to('g').magnitude
		self.cellCycleLen = kb.cellCycleLen.to('s').magnitude

		self.wildtypeBiomassReaction = kb.wildtypeBiomass['biomassFlux'].magnitude

		self.wildtypeIds = kb.wildtypeBiomass['metaboliteId']

		# Views
		self.biomassMetabolites = self.bulkMoleculesView(self.wildtypeIds)
		self.ppi = self.bulkMoleculeView("PPI[c]")
		self.ntpsdntps = self.bulkMoleculesView([
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			"DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"
			])
		self.ntps = self.bulkMoleculesView([
			"ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",
			])
		self.nmps = self.bulkMoleculesView([
			"AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"
			])

		# Attributes needed for dynamic objective

		self.nAvogadro = kb.nAvogadro.magnitude

		self.aaIdxsInWildTypeBiomass = np.array(
			[np.where(self.wildtypeIds == x)[0][0] for x in kb.aaIDs]
			)

		aaIdxsInKb = np.array([
			np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.aaIDs
			])
		self.aaMws = kb.bulkMolecules["mass"][aaIdxsInKb].magnitude

		ntpIds = ["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]",]

		self.ntpIdxsInWildTypeBiomass = np.array(
			[np.where(self.wildtypeIds == x)[0][0] for x in ntpIds]
			)

		ntpIdxsInKb = np.array([
			np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in ntpIds
			])
		self.ntpMws = kb.bulkMolecules["mass"][ntpIdxsInKb].magnitude


	def calculateRequest(self):
		self.ppi.requestAll()
		self.nmps.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		atpm = np.zeros_like(self.biomassMetabolites.counts())

		# ##### Dynamic objective #####

		requests = self.readFromListener("MetabolicDemands", "metaboliteRequests").sum(1)

		# For AAs
		relativeAArequests = normalize(requests[self.aaIdxsInWildTypeBiomass])

		if not np.any(np.isnan(relativeAArequests)):
			# print "Before: %0.10f" % (np.dot(self.biomassMws / 1000, self.wildtypeBiomassReaction))

			self.wildtypeBiomassReaction[self.aaIdxsInWildTypeBiomass] = (
				countsFromMassAndExpression(
					np.dot(
						self.aaMws / 1000,
						self.wildtypeBiomassReaction[self.aaIdxsInWildTypeBiomass]
						),
					self.aaMws,
					relativeAArequests,
					self.nAvogadro
					) *
				relativeAArequests *
				1000 / self.nAvogadro
				)

		# For NTPs
		relativeNTPrequests = normalize(requests[self.ntpIdxsInWildTypeBiomass])

		if not np.any(np.isnan(relativeNTPrequests)):
			# print "Before: %0.10f" % (np.dot(self.biomassMws / 1000, self.wildtypeBiomassReaction))

			self.wildtypeBiomassReaction[self.ntpIdxsInWildTypeBiomass] = (
				countsFromMassAndExpression(
					np.dot(
						self.ntpMws / 1000,
						self.wildtypeBiomassReaction[self.ntpIdxsInWildTypeBiomass]
						),
					self.ntpMws,
					relativeNTPrequests,
					self.nAvogadro
					) *
				relativeNTPrequests *
				1000 / self.nAvogadro
				)

			# print "After: %0.10f" % (np.dot(self.biomassMws / 1000, self.wildtypeBiomassReaction))
		# ##### End dynamic objective code #####

		deltaMetabolites = np.fmax(
			stochasticRound(
				self.randomState,
				np.round(
					(self.wildtypeBiomassReaction + atpm) * 1e-3 *
					self.nAvogadro * self.initialDryMass
					) *
				np.exp(np.log(2) / self.cellCycleLen * self.time()) *
				(np.exp(np.log(2) / self.cellCycleLen) - 1.0)
				).astype(np.int64),
			0
			)

		self.biomassMetabolites.countsInc(deltaMetabolites)

		# Fake recycling
		# if np.sum(self.ntpsdntps.counts()) < self.ppi.count():
		# 	import ipdb; ipdb.set_trace()
		# 	raise Exception, "Making fewer (d)NTPs than PPi's available"


		self.ntps.countsInc(self.nmps.counts())
		self.nmps.countsIs(0)
		self.ppi.countDec(np.sum(self.ntpsdntps.counts()))

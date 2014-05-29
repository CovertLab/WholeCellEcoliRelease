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
# import wholecell.util.flextFbaModel

from wholecell.reconstruction.fitter import normalize, countsFromMassAndExpression

class Metabolism(wholecell.processes.process.Process):
	""" Metabolism """

	_name = "Metabolism"

	# Constructor
	def __init__(self):
		# Options
		self.lpSolver = "glpk"
		self.realMax = 1e6

		# References to states
		self.metabolism = None
		self.mass = None
		self.mc = None

		# Partitions
		self.bulkMoleculesPartition = None

		# Constants
		self.initialDryMass = None
		self.cellCycleLen = None

		self.objective = None							# FBA LP objective (max growth)
		self.sMat = None								# Stoichiometry matrix [met x rxn]
		self.eMat = None								# Enzyme catalysis matrix
		self.bounds = None								# Flux bounds data
		self.dConc_dt = None							# dConc/dt (FBA LP RHS)
		self.metConstTypes = None						# Metabolite constraint types
		self.rxnVarTypes = None							# Rxn variable types
		self.rxnNewFlux = None
		self.rxnRecycFlux = None

		self.metIds = None								# Metabolite ids
		self.metNames = None							# Metabolite names
		self.metIdx = None								# Metabolite indices

		self.rxnIds = None								# Reaction ids
		self.rxnNames = None							# Reaction names
		self.rxnIdx = None								# Reaction indices

		self.biomassMetabolites = None							# Core biomass function
		self.wildtypeIds = None						# Core biomass ids

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

		### A little hacky here. John close your eyes.
		# TODO: Clean up

		# TODO: Include selenocysteine
		self.aaIds = [aaId for aaId in kb.aaIDs if aaId != "SEC-L[c]"]

		self.aaIdxsInWildTypeBiomass = np.array(
			[np.where(self.wildtypeIds == x)[0][0] for x in self.aaIds]
			)

		self.bulkMoleculesState = sim.states["BulkMolecules"]
		self.aaIdxsInContainer = self.bulkMoleculesState.container._namesToIndexes(self.aaIds)

		aaIdxsInKb = np.array([
			np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in self.aaIds
			])
		self.aaMws = kb.bulkMolecules["mass"][aaIdxsInKb].magnitude

		bulkMoleculesIdxs = np.array([
			np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in self.wildtypeIds
			])
		self.biomassMws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].magnitude # TOKB
		self.nAvogadro = kb.nAvogadro.magnitude

	def calculateRequest(self):
		self.ppi.requestAll()
		self.nmps.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		atpm = np.zeros_like(self.biomassMetabolites.counts())

		# ##### Dynamic objective #####
		relativeAArequests = normalize(
			self.bulkMoleculesState._countsRequested[self.aaIdxsInContainer].sum(axis = 1)
			)
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

			# print "After: %0.10f" % (np.dot(self.biomassMws / 1000, self.wildtypeBiomassReaction))
		# ##### End dynamic objective code #####

		deltaMetabolites = np.fmax(
			self.randStream.stochasticRound(
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


	def calcGrowthRate(self, bounds):
		growth = 1.0 / self.cellCycleLen
		fluxes = None
		exchangeRates = None

		return growth, fluxes, exchangeRates

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


	def calculateRequest(self):
		self.ppi.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		atpm = np.zeros_like(self.biomassMetabolites.counts())

		# noise = self.randStream.multivariate_normal(
		# 	np.zeros_like(self.wildtypeBiomassReaction),
		# 	np.diag(self.wildtypeBiomassReaction / 1000.)
		# 	)

		# deltaMetabolites = np.fmax(
		# 	self.randStream.stochasticRound(
		# 		np.round((self.wildtypeBiomassReaction + atpm + noise) * 1e-3
		# 			* self.nAvogadro * self.initialDryMass)
		# 		* np.exp(np.log(2) / self.cellCycleLen * self.time())
		# 		* (np.exp(np.log(2) / self.cellCycleLen) - 1.0)
		# 		).astype(np.int64),
		# 	0
		# 	)

		# No noise in production
		deltaMetabolites = np.fmax(
			self.randStream.stochasticRound(
				np.round(
					(self.wildtypeBiomassReaction) * 1e-3 *
					self.nAvogadro * self.initialDryMass
					) *
				np.exp(np.log(2) / self.cellCycleLen * self.time()) *
				(np.exp(np.log(2) / self.cellCycleLen) - 1.0)
				).astype(np.int64),
			0
			)

		self.biomassMetabolites.countsInc(deltaMetabolites)

		# Fake recycling
		if np.sum(self.ntpsdntps.counts()) < self.ppi.count():
			raise Exception, "Making fewer (d)NTPs than PPi's available"
		self.ppi.countDec(np.sum(self.ntpsdntps.counts()))


	def calcGrowthRate(self, bounds):
		growth = 1.0 / self.cellCycleLen
		fluxes = None
		exchangeRates = None

		return growth, fluxes, exchangeRates

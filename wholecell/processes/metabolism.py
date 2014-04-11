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

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "Metabolism",
			"name": "Metabolism",
			"options": ["lpSolver", "realMax"]
		}

		# Options
		self.lpSolver = "glpk"
		self.realMax = 1e6

		# References to states
		self.metabolism = None
		self.mass = None
		self.mc = None
		self.time = None

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

		self.feistCore = None							# Core biomass function
		self.feistCoreIds = None						# Core biomass ids

		super(Metabolism, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Metabolism, self).initialize(sim, kb)

		# self.mass = sim.states["Mass"]
		self.time = sim.states["Time"]

		# Load constants
		self.nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
		self.initialDryMass = kb.avgCellDryMassInit.to('g').magnitude
		self.cellCycleLen = kb.cellCycleLen.to('s').magnitude
		
		# bioIds = []
		# bioConc = []

		# for m in kb.metabolites:
		# 	if m["biomassInfo"] != {'core' : [], 'wildtype' : []}:
		# 		# TODO: Might need to fix, need to loop over all indicies within core biomass. Fix in general this is ugly.
		# 		bioIds.append("%s[%s]" % (m["id"], m["biomassInfo"]["core"][0]["location"]))
		# 		bioConc.append(m["biomassInfo"]["core"][0]["mmol/gDCW"]) # Number of molecules to produce each time step

		# bioIds, bioConc = (list(x) for x in zip(*sorted(zip(bioIds, bioConc))))
		# bioConc = np.array(bioConc)

		# self.bioProd = np.array([x if x > 0 else 0 for x in bioConc])

		self.feistCoreBiomassReaction = kb.coreBiomass['biomassFlux'].magnitude

		self.feistCoreIds = kb.coreBiomass['metaboliteId']

		# Views
		# self.biomass = self.bulkMoleculesView(bioIds)
		# self.atpHydrolysis = self.bulkMoleculesView(["ATP[c]", "H2O[c]", "ADP[c]", "PI[c]", "H[c]"])
		# self.ntps = self.bulkMoleculesView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
		# self.h2o = self.bulkMoleculeView('H2O[c]')
		self.feistCore = self.bulkMoleculesView(self.feistCoreIds)



	def calculateRequest(self):
		# self.biomass.requestAll()

		# self.ntps.requestIs(0)
		# self.h2o.requestIs(0)

		pass


	# Calculate temporal evolution
	def evolveState(self):
		atpm = np.zeros_like(self.feistCore.counts())

		noise = self.randStream.multivariate_normal(
			np.zeros_like(self.feistCoreBiomassReaction),
			np.diag(self.feistCoreBiomassReaction / 1000.)
			)

		deltaMetabolites = np.fmax(
			self.randStream.stochasticRound(
				np.round((self.feistCoreBiomassReaction + atpm + noise) * 1e-3
					* self.nAvogadro * self.initialDryMass)
				* np.exp(np.log(2) / self.cellCycleLen * self.time.value)
				* (np.exp(np.log(2) / self.cellCycleLen) - 1.0)
				).astype(np.int64),
			0
			)

		self.feistCore.countsInc(deltaMetabolites)


	def calcGrowthRate(self, bounds):
		growth = 1.0 / self.cellCycleLen
		fluxes = None
		exchangeRates = None

		return growth, fluxes, exchangeRates

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

		self.unaccountedEnergyConsumption = 6.2750e7	# ATP / cell cycle # TOKB

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
		self.nAvogadro = kb.constants['nAvogadro']['value']
		self.initialDryMass = kb.parameters['avgInitCellMass']['value'] * 10e-13 # g
		self.cellCycleLen = kb.parameters['cellCycleLen']['value']
		
		bioIds = []
		bioConc = []

		for m in kb.metabolites:
			if m["biomassInfo"] != {'core' : [], 'wildtype' : []}:
				# TODO: Might need to fix, need to loop over all indicies within core biomass. Fix in general this is ugly.
				bioIds.append("%s[%s]" % (m["id"], m["biomassInfo"]["core"][0]["location"]))
				bioConc.append(m["biomassInfo"]["core"][0]["mmol/gDCW"]) # Number of molecules to produce each time step

		bioIds, bioConc = (list(x) for x in zip(*sorted(zip(bioIds, bioConc))))
		bioConc = np.array(bioConc)

		self.bioProd = np.array([x if x > 0 else 0 for x in bioConc])

		self.feistCoreBiomassReaction = np.array([ # TODO: This needs to go in the KB
			0.513689, 0.295792, 0.241055, 0.241055, 0.091580, 0.263160, 0.263160, 0.612638, 0.094738, 0.290529,
			0.450531, 0.343161, 0.153686, 0.185265, 0.221055, 0.215792, 0.253687, 0.056843, 0.137896, 0.423162,
			0.026166, 0.027017, 0.027017, 0.026166, 0.133508, 0.215096, 0.144104, 0.174831, 0.013894, 0.019456,
			0.063814, 0.075214, 0.177645, 0.011843, 0.007895, 0.004737, 0.007106, 0.007106, 0.003158, 0.003158,
			0.003158, 0.003158, 0.003158, 0.004737, 0.003948, 0.003948, 0.000576, 0.001831, 0.000447, 0.000223,
			0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000055, 0.000223, 0.000223,
			0.000223		# mmol/gDCW (supp info 3, "biomass_core", column G)
			]) # TOKB

		self.feistCoreIds = [
			"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLN-L[c]", "GLU-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",
			"LEU-L[c]", "LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
			"DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]", "CTP[c]", "GTP[c]", "UTP[c]", "ATP[c]", "MUREIN5PX4P[p]", "KDO2LIPID4[o]",
			"PE160[c]", "PE161[c]", "K[c]", "NH4[c]", "MG2[c]", "CA2[c]", "FE2[c]", "FE3[c]", "CU2[c]", "MN2[c]",
			"MOBD[c]", "COBALT2[c]", "ZN2[c]", "CL[c]", "SO4[c]", "PI[c]", "COA[c]", "NAD[c]", "NADP[c]", "FAD[c]",
			"THF[c]", "MLTHF[c]", "10FTHF[c]", "THMPP[c]", "PYDX5P[c]", "PHEME[c]", "SHEME[c]", "UDCPDP[c]", "AMET[c]", "2OHPH[c]",
			"RIBFLV[c]"
			]

		# Views
		self.biomass = self.bulkMoleculesView(bioIds)
		self.atpHydrolysis = self.bulkMoleculesView(["ATP[c]", "H2O[c]", "ADP[c]", "PI[c]", "H[c]"])
		self.ntps = self.bulkMoleculesView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
		self.h2o = self.bulkMoleculeView('H2O[c]')
		self.feistCore = self.bulkMoleculesView(self.feistCoreIds)


	def calculateRequest(self):
		self.biomass.requestAll()

		self.ntps.requestIs(0)
		self.h2o.requestIs(0)


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

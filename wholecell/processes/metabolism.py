#!/usr/bin/env python

"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

import numpy

import wholecell.processes.process
# import wholecell.util.flextFbaModel
from wholecell.utils.constants import Constants

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
		self.mcPartition = None

		# Constants
		self.avgCellInitMass = 13.1						# fg
		self.cellCycleLen = 1.0 * 3600					# s

		self.unaccountedEnergyConsumption = 6.2750e7	# ATP / cell cycle

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
		mc = sim.states["MoleculeCounts"]
		self.time = sim.states["Time"]

		bioIds = []
		bioConc = []

		for m in kb.metabolites:
			if m["biomassConc"] != 0:
				bioIds.append("%s:%s[%s]" % (m["id"], "mature", m["biomassLoc"]))
				bioConc.append(m["biomassConc"]) # Number of molecules to produce each time step

		bioIds, bioConc = (list(x) for x in zip(*sorted(zip(bioIds, bioConc))))
		bioConc = numpy.array(bioConc)

		self.mcPartition = mc.setPartition(self, bioIds)
		self.bioProd = numpy.array([x if x > 0 else 0 for x in bioConc])

		self.mcPartition.atpHydrolysis = self.mcPartition.countsBulkViewNew(
			["ATP[c]", "H2O[c]", "ADP[c]", "PI[c]", "H[c]"])

		self.mcPartition.ntps = self.mcPartition.countsBulkViewNew(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
		self.mcPartition.h2oMol = self.mcPartition.molecule('H2O[c]')

		self.feistCore = numpy.array([ # TODO: This needs to go in the KB
			0.513689, 0.295792, 0.241055, 0.241055, 0.091580, 0.263160, 0.263160, 0.612638, 0.094738, 0.290529,
			0.450531, 0.343161, 0.153686, 0.185265, 0.221055, 0.215792, 0.253687, 0.056843, 0.137896, 0.423162,
			0.026166, 0.027017, 0.027017, 0.026166, 0.133508, 0.215096, 0.144104, 0.174831, 0.013894, 0.019456,
			0.063814, 0.075214, 0.177645, 0.011843, 0.007895, 0.004737, 0.007106, 0.007106, 0.003158, 0.003158,
			0.003158, 0.003158, 0.003158, 0.004737, 0.003948, 0.003948, 0.000576, 0.001831, 0.000447, 0.000223,
			0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000223, 0.000055, 0.000223, 0.000223,
			0.000223		# mmol/gDCW (supp info 3, "biomass_core", column G)
			])

		self.feistCoreIds = [
			"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLN-L[c]", "GLU-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",
			"LEU-L[c]", "LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
			"DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]", "CTP[c]", "GTP[c]", "UTP[c]", "ATP[c]", "MUREIN5PX4P[p]", "KDO2LIPID4[o]",
			"PE160[c]", "PE161[c]", "K[c]", "NH4[c]", "MG2[c]", "CA2[c]", "FE2[c]", "FE3[c]", "CU2[c]", "MN2[c]",
			"MOBD[c]", "COBALT2[c]", "ZN2[c]", "CL[c]", "SO4[c]", "PI[c]", "COA[c]", "NAD[c]", "NADP[c]", "FAD[c]",
			"THF[c]", "MLTHF[c]", "10FTHF[c]", "THMPP[c]", "PYDX5P[c]", "PHEME[c]", "SHEME[c]", "UDCPDP[c]", "AMET[c]", "2OHPH[c]",
			"RIBFLV[c]"
			]

		self.mcPartition.feistCore = self.mcPartition.countsBulkViewNew(self.feistCoreIds)

		self.initialDryMass = 2.8e-13 / 1.36 # grams

	# Calculate needed metabolites
	def requestMoleculeCounts(self):
		self.mcPartition.countsBulkIs(1)

		self.mcPartition.ntps.countsBulkIs(0)
		self.mcPartition.h2oMol.countBulkIs(0)

	# # Calculate needed proteins
	# def calcReqEnzyme(self):
	# 	return numpy.ones(self.enzyme.fullCounts.shape)

	# Calculate temporal evolution
	def evolveState(self):
		# NOTE: I've deleted a bunch of commented-out code from here.  Get it 
		# from an old commit if you really need it. - JM

		atpm = numpy.zeros_like(self.mcPartition.feistCore.countsBulk()) # TODO: determine what this means

		noise = self.randStream.multivariate_normal(numpy.zeros_like(self.feistCore), numpy.diag(self.feistCore / 1000.))

		deltaMetabolites = (
			numpy.round((self.feistCore + atpm + noise) * 1e-3
				* Constants.nAvogadro * self.initialDryMass)
			* numpy.exp(numpy.log(2) / self.cellCycleLen * self.time.value)
			* (numpy.exp(numpy.log(2) / self.cellCycleLen) - 1.0)
			)

		self.mcPartition.feistCore.countsBulkIs(
			numpy.fmax(
				0,
				self.randStream.stochasticRound(self.mcPartition.feistCore.countsBulk() + deltaMetabolites)
				)
			)


	def calcGrowthRate(self, bounds):
		growth = 1.0 / self.cellCycleLen
		fluxes = None
		exchangeRates = None

		return growth, fluxes, exchangeRates

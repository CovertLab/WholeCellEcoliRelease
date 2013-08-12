#!/usr/bin/env python

"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

import numpy

import wholecell.sim.process.Process

class Metabolism(wholecell.sim.process.Process.Process):
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
		self.metabolite = None
		self.enzyme = None
		self.mass = None
		self.mc = None
		self.time = None

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

		super(Metabolism, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Metabolism, self).initialize(sim, kb)

		self.mass = sim.getState("Mass")
		self.mc = sim.getState("MoleculeCounts")
		self.time = sim.getState("Time")

		bioIds = []
		bioConc = []

		for m in kb.metabolites:
			if m["biomassConc"] != 0:
				bioIds.append("%s:%s[%s]" % (m["id"], "mature", m["biomassLoc"]))
				bioConc.append(m["biomassConc"]) # Number of molecules to produce each time step

		bioIds, bioConc = (list(x) for x in zip(*sorted(zip(bioIds, bioConc))))
		bioConc = numpy.array(bioConc)

		self.metabolite = sim.getState("MoleculeCounts").addPartition(self, bioIds, self.calcReqMetabolites)
		self.metabolite.idx["atpHydrolysis"] = self.metabolite.getIndex(["ATP[c]", "H2O[c]", "ADP[c]", "PI[c]", "H[c]"])[0]
		self.metabolite.idx["ntps"] = self.metabolite.getIndex(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])[0]
		self.metabolite.idx["h2o"]  = self.metabolite.getIndex("H2O[c]")[0]
		self.bioProd = numpy.array([x if x > 0 else 0 for x in bioConc])

		self.metabolite.idx["FeistCoreRows"], self.metabolite.idx["FeistCoreCols"] = self.metabolite.getIndex([
			"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLN-L[c]", "GLU-L[c]", "GLY[c]", "HIS-L[c]", "ILE-L[c]",
			"LEU-L[c]", "LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]",
			"DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]", "CTP[c]", "GTP[c]", "UTP[c]", "ATP[c]", "MUREIN5PX4P[p]", "KDO2LIPID4[o]",
			"PE160[c]", "PE161[c]", "K[c]", "NH4[c]", "MG2[c]", "CA2[c]", "FE2[c]", "FE3[c]", "CU2[c]", "MN2[c]",
			"MOBD[c]", "COBALT2[c]", "ZN2[c]", "CL[c]", "SO4[c]", "PI[c]", "COA[c]", "NAD[c]", "NADP[c]", "FAD[c]",
			"THF[c]", "MLTHF[c]", "10FTHF[c]", "THMPP[c]", "PYDX5P[c]", "PHEME[c]", "SHEME[c]", "UDCPDP[c]", "AMET[c]", "2OHPH[c]",
			"RIBFLV[c]"
			])[1:]

	# Calculate needed metabolites
	def calcReqMetabolites(self):
		val = numpy.ones(self.metabolite.fullCounts.shape)
		val[self.metabolite.idx["ntps"]] = 0
		val[self.metabolite.idx["h2o"]] = 0
		return val

	# Calculate needed proteins
	def calcReqEnzyme(self):
		return numpy.ones(self.enzyme.fullCounts.shape)

	# Calculate temporal evolution
	def evolveState(self):
		from wholecell.util.Constants import Constants
		dm = numpy.zeros(self.metabolite.counts.shape)
		atpm = numpy.zeros(self.metabolite.counts.shape)
		atpm[self.metabolite.idx["ntps"][0]] = 12.4
		noise = numpy.random.multivariate_normal(numpy.zeros(self.mc.vals["FeistCore"].size), numpy.diag(self.mc.vals["FeistCore"] / 100.))
		dm[self.metabolite.idx["FeistCoreRows"]] = numpy.round((self.mc.vals["FeistCore"] + atpm[self.metabolite.idx["FeistCoreRows"]] + noise ) * 1e-3 * Constants.nAvogadro * self.mc.initialDryMass) * numpy.exp(numpy.log(2) / self.cellCycleLen * self.time.value) * (numpy.exp(numpy.log(2) / self.cellCycleLen) - 1.0)
		# import ipdb
		# ipdb.set_trace()
		print "NTP production: %s" % str(dm[self.metabolite.idx["ntps"]])
		self.metabolite.counts = self.randStream.stochasticRound(
			self.metabolite.counts + dm
			)

		# Unaccounted energy consumption
		# self.metabolite.counts[self.metabolite.idx["atpHydrolysis"]] = \
		# 	+ self.metabolite.counts[self.metabolite.idx["atpHydrolysis"]] \
		# 	+ numpy.array([-1.0, -1.0, 1.0, 1.0, 1.0]) * self.randStream.stochasticRound(self.unaccountedEnergyConsumption * growth_cellPerSec * self.timeStepSec)

		# Make copy numbers positive
		self.metabolite.counts = numpy.maximum(0, self.metabolite.counts)

	def calcGrowthRate(self, bounds):
		growth = 1.0 / self.cellCycleLen
		fluxes = None
		exchangeRates = None

		return growth, fluxes, exchangeRates
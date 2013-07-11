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

		bioIds = []
		bioConc = []

		for m in kb.metabolites:
			if m["biomassConc"] > 0:
				bioIds.append("%s:%s[%s]" % (m["id"], "mature", m["biomassLoc"]))
				bioConc.append(self.timeStepSec * m["biomassConc"] / self.cellCycleLen) # Number of molecules to produce each time step

		bioIds, bioConc = (list(x) for x in zip(*sorted(zip(bioIds, bioConc))))
		bioConc = numpy.array(bioConc)

		self.metabolite = sim.getState("MoleculeCounts").addPartition(self, bioIds, self.calcReqMetabolites)
		self.bioProd = bioProd

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
		self.metabolite.counts = self.randStream.stochasticRound(
			self.metabolite.counts + self.bioProd
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
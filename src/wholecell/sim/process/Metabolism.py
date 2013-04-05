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
		self.cellCycleLen = 9 * 3600					# s

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

		# List of metabolites, enzymes
		molIds = []
		enzIds = []
		for r in kb.reactions:
			for s in r["stoichiometry"]:
				molIds.append("%s:%s[%s]" % (s["molecule"], s["form"], s["compartment"]))
			e = r["enzyme"]
			if type(e) == dict:
				enzIds.append("%s:%s[%s]" % (e["id"], e["form"], e["compartment"]))
		molIds = list(set(molIds))
		enzIds = list(set(enzIds))

		# Partitions
		self.metabolism = sim.getState("Metabolism").addPartition(self)
		self.metabolite = sim.getState("MoleculeCounts").addPartition(self, molIds, self.calcReqMetabolites)
		self.enzyme = sim.getState("MoleculeCounts").addPartition(self, enzIds, self.calcReqEnzyme)
		self.mass = sim.getState("Mass").addPartition(self)

		self.metabolite.idx["atpHydrolysis"] = self.metabolite.getIndex(["ATP[c]", "H2O[c]", "ADP[c]", "PI[c]", "H[c]"])[0]
		self.metabolite.idx["ntps"] = self.metabolite.getIndex(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])[0]
		self.metabolite.idx["ndps"] = self.metabolite.getIndex(["ADP[c]", "CDP[c]", "GDP[c]", "UDP[c]"])[0]
		self.metabolite.idx["nmps"] = self.metabolite.getIndex(["AMP[c]", "CMP[c]", "GMP[c]", "UMP[c]"])[0]
		self.metabolite.idx["ppi"]  = self.metabolite.getIndex("PPI[c]")[0]
		self.metabolite.idx["pi"]   = self.metabolite.getIndex("PI[c]")[0]
		self.metabolite.idx["h2o"]  = self.metabolite.getIndex("H2O[c]")[0]
		self.metabolite.idx["h"]    = self.metabolite.getIndex("H[c]")[0]

		# Indices
		nExchangeConstraints = 7

		nMet = len(molIds) + 1 + nExchangeConstraints
		self.metIds = self.metabolite.ids + ["biomass"]
		self.metNames = self.metabolite.ids + ["biomass"]
		self.metIdx = {}
		self.metIdx["real"] = numpy.arange(len(molIds))				# Real metabolites
		self.metIdx["biomass"] = numpy.array([self.metIdx["real"][-1]]) + 1
		self.metIdx["exchangeConstraints"] = numpy.array(self.metIdx["biomass"][-1] + numpy.arange(nExchangeConstraints) + 1)

		nRxn = len(kb.reactions) + len(molIds) + 2
		mc = sim.getState("MoleculeCounts")
		cIdxs = numpy.unravel_index(self.metabolite.mapping, (len(mc.ids), len(mc.compartments)))[1]
		self.rxnIds = \
			self.metabolism.reactionIds + \
			["ex_" + x[0] + "_" + x[1] for x in zip(self.metabolite.ids, [mc.compartments[x]["id"] for x in cIdxs])] + \
			["growth"] + \
			["ex_biomass"]
		self.rxnNames = \
			self.metabolism.reactionNames + \
			["ex_" + x[0] + "_" + x[1] for x in zip(self.metabolite.names, [mc.compartments[x]["name"] for x in cIdxs])] + \
			["growth"] + \
			["biomass exchange"]
		self.rxnIdx = {}
		self.rxnIdx["real"] = numpy.arange(len(kb.reactions))
		self.rxnIdx["exchange"] = self.rxnIdx["real"][-1] + numpy.arange(len(molIds)) + 1
		self.rxnIdx["growth"] = numpy.array([self.rxnIdx["exchange"][-1]]) + 1
		self.rxnIdx["biomassExchange"] = numpy.array(self.rxnIdx["growth"][-1] + 1)

		mc = sim.getState("MoleculeCounts")
		cIdxs = numpy.array(mc.getIndex(molIds)[2])
		iExtracellular = [x["id"] for x in mc.compartments].index("e")
		self.rxnIdx["internalExchange"] = self.rxnIdx["exchange"][cIdxs != iExtracellular]
		self.rxnIdx["externalExchange"] = self.rxnIdx["exchange"][cIdxs == iExtracellular]

		nEnz = len(enzIds)

		# Stoichiometry matrix, enzymes, kinetics
		self.sMat = numpy.zeros([nMet, nRxn])
		self.sMat[self.metIdx["real"], self.rxnIdx["exchange"]] = numpy.ones(len(molIds))
		self.sMat[self.metIdx["biomass"], self.rxnIdx["growth"]] = 1.0
		self.sMat[self.metIdx["biomass"], self.rxnIdx["biomassExchange"]] = -1.0

		self.eMat = numpy.zeros([nRxn, nEnz])
		self.bounds = {
			"kinetic": {"lo": -numpy.ones(nRxn) * numpy.Inf, "up": numpy.ones(nRxn) * numpy.Inf},
			"thermodynamic": {"lo": -numpy.ones(nRxn) * numpy.Inf, "up": numpy.ones(nRxn) * numpy.Inf},
			"exchange": {"lo": -numpy.ones(nRxn) * numpy.Inf, "up": numpy.ones(nRxn) * numpy.Inf}
		}
		self.bounds["thermodynamic"]["lo"][self.rxnIdx["growth"]] = 0
		self.bounds["thermodynamic"]["lo"][self.rxnIdx["biomassExchange"]] = 0

		tmpSMat = []
		tmpEMat = []
		for iReaction in xrange(len(kb.reactions)):
			r = kb.reactions[iReaction]
			rIdx = self.rxnIdx["real"][iReaction]

			# stoichiometry
			for s in r["stoichiometry"]:
				tmpSMat.append(["%s:%s[%s]" % (s["molecule"], s["form"], s["compartment"]), rIdx, s["coeff"]])

			# enzyme
			if type(r["enzyme"]) == dict:
				# Catalysis
				e = r["enzyme"]
				tmpEMat.append(["%s:%s[%s]" % (e["id"], e["form"], e["compartment"]), rIdx])

				# Kinetics
				if not numpy.isnan(e["kCatRev"]):
					self.bounds["kinetic"]["lo"][rIdx] = -e["kCatRev"]
				if not numpy.isnan(e["kCatFor"]):
					self.bounds["kinetic"]["up"][rIdx] = e["kCatFor"]

			# thermodynamics
			if r["dir"] == 1:
				self.bounds["thermodynamic"]["lo"][rIdx] = 0
			elif r["dir"] == -1:
				self.bounds["thermodynamic"]["up"][rIdx] = 0
		mIdx = self.metIdx["real"][numpy.array(self.metabolite.getIndex([x[0] for x in tmpSMat])[0])]
		self.sMat[mIdx, numpy.array([x[1] for x in tmpSMat])] = numpy.array([x[2] for x in tmpSMat])

		eIdx = numpy.array(self.enzyme.getIndex([x[0] for x in tmpEMat])[0])
		self.eMat[[x[1] for x in tmpEMat], eIdx] = 1.0

		# exchange
		metIds = [x["id"] + ":mature[e]" for x in kb.metabolites if x["id"] + ":mature[e]" in molIds]
		metExs = numpy.array([x["maxExchangeRate"] for x in kb.metabolites if x["id"] in molIds])

		metIdxs = numpy.array(self.metabolite.getIndex(metIds)[0])
		self.bounds["exchange"]["lo"][self.rxnIdx["exchange"][metIdxs]] = -metExs
		self.bounds["exchange"]["up"][self.rxnIdx["exchange"][metIdxs]] =  metExs

		# exchange constraints
		self.sMat[self.metIdx["exchangeConstraints"][0], self.rxnIdx["exchange"][numpy.array(self.metabolite.getIndex(["ATP[c]", "ADP[c]", "AMP[c]"])[0])]] = 1.0
		self.sMat[self.metIdx["exchangeConstraints"][1], self.rxnIdx["exchange"][numpy.array(self.metabolite.getIndex(["CTP[c]", "CDP[c]", "CMP[c]"])[0])]] = 1.0
		self.sMat[self.metIdx["exchangeConstraints"][2], self.rxnIdx["exchange"][numpy.array(self.metabolite.getIndex(["GTP[c]", "GDP[c]", "GMP[c]"])[0])]] = 1.0
		self.sMat[self.metIdx["exchangeConstraints"][3], self.rxnIdx["exchange"][numpy.array(self.metabolite.getIndex(["UTP[c]", "UDP[c]", "UMP[c]"])[0])]] = 1.0
		self.sMat[self.metIdx["exchangeConstraints"][4], self.rxnIdx["exchange"][numpy.array(self.metabolite.getIndex(["FTHF10[c]", "THF[c]"])[0])]] = 1.0
		self.sMat[self.metIdx["exchangeConstraints"][5], self.rxnIdx["exchange"][numpy.array(self.metabolite.getIndex(["FTHF10[c]", "FOR[c]", "FMET[c]"])[0])]] = 1.0
		self.sMat[self.metIdx["exchangeConstraints"][6], self.rxnIdx["exchange"][numpy.array(self.metabolite.getIndex(["MET[c]", "FMET[c]"])[0])]] = 1.0

		# objective
		objMets = [x for x in kb.metabolites if x["metabolismNewFlux"] != 0 or x["metabolismRecyclingFlux"] != 0]

		metComps = ["m" if x["hydrophobic"] == True else "c" for x in objMets]

		realMetIds = [x[0] + "[" + x[1] + "]" for x in zip([x["id"] for x in objMets], metComps)]
		realMetIdxs = numpy.array(self.metabolite.getIndex(realMetIds)[0])

		self.rxnNewFlux = numpy.zeros(nMet)
		self.rxnRecycFlux = numpy.zeros(nMet)
		self.rxnNewFlux[self.metIdx["real"][realMetIdxs]] = numpy.array([x[0]["metabolismNewFlux"] for x in zip(objMets, realMetIdxs) if x[1] != 0])
		self.rxnRecycFlux[self.metIdx["real"][realMetIdxs]] = numpy.array([x[0]["metabolismRecyclingFlux"] for x in zip(objMets, realMetIdxs) if x[1] != 0])

		self.rxnIdx["internalNoRecycExchange"] = numpy.intersect1d(self.rxnIdx["exchange"][self.rxnRecycFlux[self.rxnIdx["real"]] == 0], self.rxnIdx["internalExchange"])
		self.rxnIdx["internalRecycExchange"] = numpy.intersect1d(self.rxnIdx["exchange"][self.rxnRecycFlux[self.rxnIdx["real"]] < 0], self.rxnIdx["internalExchange"])

		self.sMat[self.metIdx["real"], self.rxnIdx["growth"]] = -self.rxnNewFlux[self.metIdx["real"]]

		self.objective = numpy.zeros(nRxn)
		self.objective[self.rxnIdx["growth"]] = 1e3
		self.objective[self.rxnIdx["internalRecycExchange"]] = 1.0 / numpy.sum(numpy.minimum(0, self.rxnNewFlux))

		# dConc/dt
		self.dConc_dt = numpy.zeros(nMet)

		# constraint, variable types
		self.metConstTypes = ["S"] * nMet
		self.rxnVarTypes = ["C"] * nRxn

		# more indices
		self.rxnIdx["catalyzed"] = numpy.where(numpy.any(self.eMat, axis = 1))[0]

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
		# Calculate flux bounds
		bounds = self.calcFluxBounds[self.metabolite.counts, self.enzyme.counts]

		# Calculate growth rate
		self.metabolism.growth, self.metabolism.fluxes, exRates = self.calcGrowthRate(bounds)
		growth_cellPerSec = self.metabolism.growth / (3600.0 * self.avgCellInitMass)

		# Update metabolite copy numbers
		self.metabolite.counts = self.randStream.stochasticRound(
			+ self.metabolite.counts
			- self.sMat[self.metIdx["real"], self.rxnIdx["growth"]] * growth_cellPerSec * self.timeStepSec
			+ exRates * self.timeStepSec
			)

		# Unaccounted energy consumption
		self.metabolite.counts[self.metabolite.idx["atpHydrolysis"]] = \
			+ self.metabolite.counts[self.metabolite.idx["atpHydrolysis"]] \
			+ numpy.array([-1.0, -1.0, 1.0, 1.0, 1.0]) * self.randStream.stochasticRound(self.unaccountedEnergyConsumption * growth_cellPerSec * self.timeStepSec)

		# Make copy numbers positive
		self.metabolite.counts = numpy.max(0, self.metabolite.counts)

	def calcGrowthRate(self, bounds):
		import wholecell.util.linearProgramming as linearProgramming

		# Cap bounds
		bounds["lo"] = numpy.mininum(0, numpy.maximum(bounds["lo"] * self.timeStepSec, -self.realMax))
		bounds["up"] = numpy.maximum(0, numpy.minimum(bounds["up"] * self.timeStepSec,  self.realMax))

		# Flux-balance analysis
		x, sol = linearProgramming.linearProgramming(
			"maximize", self.objective, self.sMat, self.dConc_dt,
			bounds["lo"], bounds["up"], self.metConstTypes, self.rxnVarTypes,
			None	# Will use glpk
			)

		if sol["status"] != "optimal":
			raise Exception, "Linear programming error. Status: %s" % sol["status"]

		# extract growth, real fluxes
		x /= self.timeStepSec
		growth = x[self.rxnIdx["growth"]] * 3600 * self.avgCellInitMass
		fluxes = x[self.rxnIdx["real"]]
		exchangeRates = x[self.rxnIdx["exchange"]]

		return growth, fluxes, exchangeRates

	def calcFluxBounds(self, metCnts, enzCnts,
		applyThermoBounds = True, applyKineticBounds = True,
		applyExchangeBounds = True, applyMetAvailBounds = True):

		# Initialize
		nRxn = len(self.rxnIds)
		lo = -numpy.ones(nRxn) * numpy.Inf
		up =  numpy.ones(nRxn) * numpy.Inf

		# Thermodynamics
		if applyThermoBounds:
			lo = numpy.maximum(lo, self.bounds["thermodynamic"]["lo"])
			up = numpy.minimum(up, self.bounds["thermodynamic"]["up"])

		# Kinetics
		if applyKineticBounds:
			rxnEnzCnts = numpy.dot(self.eMat[self.rxnIdx["catalyzed"], :], enzCnts)

			kLo = rxnEnzCnts * self.bounds["kinetic"]["lo"][self.rxnIdx["catalyzed"]]
			kUp = rxnEnzCnts * self.bounds["kinetic"]["up"][self.rxnIdx["catalyzed"]]
			kLo[rxnEnzCnts == 0] = 0
			kUp[rxnEnzCnts == 0] = 0

			lo[self.rxnIdx["catalyzed"]] = numpy.maximum(lo[self.rxnIdx["catalyzed"]], kLo)
			up[self.rxnIdx["catalyzed"]] = numpy.minimum(up[self.rxnIdx["catalyzed"]], kUp)

		# Exchange
		if applyExchangeBounds:
			lo[self.rxnIdx["internalNoRecycExchange"]] = 0
			up[self.rxnIdx["internalNoRecycExchange"]] = 0

			lo[self.rxnIdx["externalExchange"]] = numpy.maximum(lo[self.rxnIdx["externalExchange"]],
				self.bounds["exchange"]["lo"][self.rxnIdx["externalExchange"]] * 6.022e23 * 1e-3 * 3600 * numpy.sum(self.mass.cellDry) * 1e-15)
			up[self.rxnIdx["externalExchange"]] = numpy.minimum(up[self.rxnIdx["externalExchange"]],
				self.bounds["exchange"]["up"][self.rxnIdx["externalExchange"]] * 6.022e23 * 1e-3 * 3600 * numpy.sum(self.mass.cellDry) * 1e-15)

		# Metabolite availability
		if applyMetAvailBounds:
			lo[self.rxnIdx["exchange"]] = numpy.maximum(lo[self.rxnIdx["exchange"]], -metCnts / self.timeStepSec)

		# Protein # TODO

		return {"lo": lo, "up": up}
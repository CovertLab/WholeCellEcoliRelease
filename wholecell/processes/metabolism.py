#!/usr/bin/env python

"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

TODO:
- move over to flexFBA
- implement metabolite pools
- enzyme-limited reactions (& fit enzyme expression)
- option to call a reduced form of metabolism (assume optimal)

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/2/2013
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
import wholecell.utils.flex_t_fba_model
import wholecell.utils.d_fba_model

from wholecell.reconstruction.fitter import normalize, countsFromMassAndExpression
from wholecell.utils.random import stochasticRound
from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM

# from wholecell.utils.modular_fba import FluxBalanceAnalysis
from wholcell.utils.modular_fba import _setupFeist as setupFeist # TODO: move setup to KB/process

GAIN = 10
ASSUME_OPTIMAL_GROWTH = True

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

		self.wildtypeBiomassReaction = (
			kb.wildtypeBiomass['biomassFlux'].magnitude
			)
		self.wildtypeBiomassPoolIncreases = (
			kb.wildtypeBiomassPoolIncreases["biomassFlux"].magnitude *
			self.timeStepSec
			)
		self.wildtypeBiomassReactionSS = self.wildtypeBiomassReaction.copy()

		self.wildtypeIds = kb.wildtypeBiomass['metaboliteId']

		# Views
		self.biomassMetabolites = self.bulkMoleculesView(self.wildtypeIds)
		self.ppi = self.bulkMoleculeView("PPI[c]")
		self.atpRecyclingReactants = self.bulkMoleculesView(
			["ADP[c]", "PI[c]", "H[c]"]
			)
		self.atpRecyclingProducts = self.bulkMoleculesView(
			["ATP[c]", "H2O[c]"]
			)
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
		self.aas = self.bulkMoleculesView(kb.aaIDs)
		self.dntps = self.bulkMoleculesView([
			"DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"
			])

		self.pi = self.bulkMoleculeView("PI[c]")
		self.h = self.bulkMoleculeView("H[c]")
		self.h2o = self.bulkMoleculeView("H2O[c]")

		# Attributes needed for dynamic objective

		self.nAvogadro = kb.nAvogadro.magnitude

		self.aaIdxsInWildTypeBiomass = np.array(
			[np.where(self.wildtypeIds == x)[0][0] for x in kb.aaIDs]
			)

		aaIdxsInKb = np.array([
			np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in kb.aaIDs
			])
		self.aaMws = kb.bulkMolecules["mass"][aaIdxsInKb].magnitude.sum(1)


		ntpIds = ["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"]

		self.ntpIdxsInWildTypeBiomass = np.array(
			[np.where(self.wildtypeIds == x)[0][0] for x in ntpIds]
			)

		ntpIdxsInKb = np.array([
			np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in ntpIds
			])
		self.ntpMws = kb.bulkMolecules["mass"][ntpIdxsInKb].magnitude.sum(1)


		dntpIds = ["DATP[c]", "DCTP[c]", "DGTP[c]", "DTTP[c]"]

		self.dntpIdxsInWildTypeBiomass = np.array(
			[np.where(self.wildtypeIds == x)[0][0] for x in dntpIds]
			)

		self.notDntpIdxsInWildTypeBiomass = np.array(
			[np.where(self.wildtypeIds == x)[0][0] for x in self.wildtypeIds if x not in dntpIds]
			)

		self.notAADntpNtpIdxsInWildTypeBiomass = np.array(
			[np.where(self.wildtypeIds == x)[0][0] for x in self.wildtypeIds if x not in dntpIds and x not in ntpIds and x not in kb.aaIDs]
			)

		dntpIdxsInKb = np.array([
			np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in dntpIds
			])

		self.dntpMws = kb.bulkMolecules["mass"][dntpIdxsInKb].magnitude.sum(1)

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)

		# Please don't delete this. It is useful for debugging.
		bulkMoleculesIdxs = np.array([
			np.where(kb.bulkMolecules["moleculeId"] == x)[0][0] for x in self.wildtypeIds
			])
		self.biomassMws = kb.bulkMolecules["mass"][bulkMoleculesIdxs].magnitude

		# biomass = [{'coeff': -0.000223, 'id': '10fthf[c]'}, {'coeff': -0.000223, 'id': '2ohph[c]'}, {'coeff': 59.810000000000002, 'id': 'adp[c]'}, {'coeff': -0.51370000000000005, 'id': 'ala-L[c]'}, {'coeff': -0.000223, 'id': 'amet[c]'}, {'coeff': -0.29580000000000001, 'id': 'arg-L[c]'}, {'coeff': -0.24110000000000001, 'id': 'asn-L[c]'}, {'coeff': -0.24110000000000001, 'id': 'asp-L[c]'}, {'coeff': -59.984000000000002, 'id': 'atp[c]'}, {'coeff': -0.0047369999999999999, 'id': 'ca2[c]'}, {'coeff': -0.0047369999999999999, 'id': 'cl[c]'}, {'coeff': -0.00057600000000000001, 'id': 'coa[c]'}, {'coeff': -0.0031580000000000002, 'id': 'cobalt2[c]'}, {'coeff': -0.13350000000000001, 'id': 'ctp[c]'}, {'coeff': -0.0031580000000000002, 'id': 'cu2[c]'}, {'coeff': -0.091579999999999995, 'id': 'cys-L[c]'}, {'coeff': -0.026169999999999999, 'id': 'datp[c]'}, {'coeff': -0.027019999999999999, 'id': 'dctp[c]'}, {'coeff': -0.027019999999999999, 'id': 'dgtp[c]'}, {'coeff': -0.026169999999999999, 'id': 'dttp[c]'}, {'coeff': -0.000223, 'id': 'fad[c]'}, {'coeff': -0.0071060000000000003, 'id': 'fe2[c]'}, {'coeff': -0.0071060000000000003, 'id': 'fe3[c]'}, {'coeff': -0.26319999999999999, 'id': 'gln-L[c]'}, {'coeff': -0.26319999999999999, 'id': 'glu-L[c]'}, {'coeff': -0.61260000000000003, 'id': 'gly[c]'}, {'coeff': -0.21510000000000001, 'id': 'gtp[c]'}, {'coeff': -54.462000000000003, 'id': 'h2o[c]'}, {'coeff': 59.810000000000002, 'id': 'h[c]'}, {'coeff': -0.094740000000000005, 'id': 'his-L[c]'}, {'coeff': -0.29049999999999998, 'id': 'ile-L[c]'}, {'coeff': -0.17760000000000001, 'id': 'k[c]'}, {'coeff': -0.019449999999999999, 'id': 'kdo2lipid4[e]'}, {'coeff': -0.45050000000000001, 'id': 'leu-L[c]'}, {'coeff': -0.34320000000000001, 'id': 'lys-L[c]'}, {'coeff': -0.1537, 'id': 'met-L[c]'}, {'coeff': -0.0078949999999999992, 'id': 'mg2[c]'}, {'coeff': -0.000223, 'id': 'mlthf[c]'}, {'coeff': -0.0031580000000000002, 'id': 'mn2[c]'}, {'coeff': -0.0031580000000000002, 'id': 'mobd[c]'}, {'coeff': -0.01389, 'id': 'murein5px4p[p]'}, {'coeff': -0.0018309999999999999, 'id': 'nad[c]'}, {'coeff': -0.00044700000000000002, 'id': 'nadp[c]'}, {'coeff': -0.011842999999999999, 'id': 'nh4[c]'}, {'coeff': -0.022329999999999999, 'id': 'pe160[c]'}, {'coeff': -0.041480000000000003, 'id': 'pe160[p]'}, {'coeff': -0.02632, 'id': 'pe161[c]'}, {'coeff': -0.048890000000000003, 'id': 'pe161[p]'}, {'coeff': -0.1759, 'id': 'phe-L[c]'}, {'coeff': -0.000223, 'id': 'pheme[c]'}, {'coeff': 59.805999999999997, 'id': 'pi[c]'}, {'coeff': 0.77390000000000003, 'id': 'ppi[c]'}, {'coeff': -0.22109999999999999, 'id': 'pro-L[c]'}, {'coeff': -0.000223, 'id': 'pydx5p[c]'}, {'coeff': -0.000223, 'id': 'ribflv[c]'}, {'coeff': -0.21579999999999999, 'id': 'ser-L[c]'}, {'coeff': -0.000223, 'id': 'sheme[c]'}, {'coeff': -0.0039480000000000001, 'id': 'so4[c]'}, {'coeff': -0.000223, 'id': 'thf[c]'}, {'coeff': -0.000223, 'id': 'thmpp[c]'}, {'coeff': -0.25369999999999998, 'id': 'thr-L[c]'}, {'coeff': -0.056840000000000002, 'id': 'trp-L[c]'}, {'coeff': -0.13789999999999999, 'id': 'tyr-L[c]'}, {'coeff': -5.5000000000000002e-05, 'id': 'udcpdp[c]'}, {'coeff': -0.14410000000000001, 'id': 'utp[c]'}, {'coeff': -0.42320000000000002, 'id': 'val-L[c]'}, {'coeff': -0.0031580000000000002, 'id': 'zn2[c]'}]
		# for x in biomass:
		# 	x["id"] = x["id"][:-2].upper() + x["id"][-2:]
		# 	if x["id"] == "KDO2LIPID4[e]":
		# 		x["id"] = "KDO2LIPID4[o]"


		# metIds = [x for x in kb.metabolismMoleculeNames]
		# import re
		# rxns = [x for x in kb.metabolismBiochemicalReactions if not re.match(".*_[0-9]$", x["id"]) or x["id"].endswith("_0") or "PFK_2" in x["id"]]
		# rxnIsIrreversible = np.array([x["dir"] for x in rxns], dtype = np.bool)
		# mediaEx = kb.metabolismMediaEx
		# # biomass = [{"id": x["metaboliteId"], "coeff": -x["biomassFlux"]} for x in kb.wildtypeBiomass.struct_array]
		# atpId = "ATP[c]"
		# # self.flexTFbaModel = wholecell.utils.flex_t_fba_model.FlexTFbaModel(metIds = metIds, rxns = rxns, mediaEx = mediaEx, biomass = biomass, atpId = "ATP[c]", params = None)
		# self.flexTFbaModel = wholecell.utils.d_fba_model.dFbaModel(metIds = metIds, rxns = rxns, mediaEx = mediaEx, biomass = biomass)
		# self.lb = wholecell.utils.flex_t_fba_model.bounds(["thermodynamic", "exchange", "bs"], self.flexTFbaModel.rxnIds(), False)
		# self.ub = wholecell.utils.flex_t_fba_model.bounds(["thermodynamic", "exchange", "bs"], self.flexTFbaModel.rxnIds(), True)

		# # Thermodynamic bounds
		# self.lb.valuesIs(self.flexTFbaModel.rxnGroup("real").idxs()[rxnIsIrreversible], "thermodynamic", 0)

		# # Biomass return is zero (for testing)
		# # self.lb.valuesIs(self.flexTFbaModel.rxnGroup("x").idxs(), "exchange", 0)

		# # Media exchange
		# self.lb.valuesIs(self.flexTFbaModel.rxnGroup("mediaEx").idxs(), "exchange", 0)
		# mediaRxnIds = [
		# 	"mediaEx_FEIST_EX_ca2(e)",
		# 	"mediaEx_FEIST_EX_cl(e)",
		# 	"mediaEx_FEIST_EX_co2(e)",
		# 	"mediaEx_FEIST_EX_cobalt2(e)",
		# 	"mediaEx_FEIST_EX_cu2(e)",
		# 	"mediaEx_FEIST_EX_fe2(e)",
		# 	"mediaEx_FEIST_EX_fe3(e)",
		# 	"mediaEx_FEIST_EX_h(e)",
		# 	"mediaEx_FEIST_EX_h2o(e)",
		# 	"mediaEx_FEIST_EX_k(e)",
		# 	"mediaEx_FEIST_EX_mg2(e)",
		# 	"mediaEx_FEIST_EX_mn2(e)",
		# 	"mediaEx_FEIST_EX_mobd(e)",
		# 	"mediaEx_FEIST_EX_na1(e)",
		# 	"mediaEx_FEIST_EX_nh4(e)",
		# 	"mediaEx_FEIST_EX_pi(e)",
		# 	"mediaEx_FEIST_EX_so4(e)",
		# 	"mediaEx_FEIST_EX_tungs(e)",
		# 	"mediaEx_FEIST_EX_zn2(e)",
		# 	# "mediaEx_FEIST_EX_cbl1(e)",
		# 	# "mediaEx_SELNP_MEDIA_EXCHANGE_HACKED",
		# 	]
		# self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(mediaRxnIds), "exchange", -np.inf)
		# self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["mediaEx_FEIST_EX_cbl1(e)"]), "exchange", -0.01)

		# # Nutrient limitations
		# self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["mediaEx_FEIST_EX_glc(e)"]), "exchange", -8.)
		# self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["mediaEx_FEIST_EX_o2(e)"]), "exchange", -18.5)
		# # self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["mediaEx_FEIST_EX_glc(e)"]), "exchange", -100.)
		# # self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["mediaEx_FEIST_EX_o2(e)"]), "exchange", -100)

		# # ATPM
		# self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_ATPM"]), "bs", 8.39)
		# self.ub.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_ATPM"]), "bs", 8.39)

		# # Reactions Feist arbitrarily set to 0
		# self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_CAT_0"]), "bs", 0)
		# self.ub.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_CAT_0"]), "bs", 0)
		# # self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_CAT_1"]), "bs", 0)
		# # self.ub.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_CAT_1"]), "bs", 0)
		# self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_SPODM_0"]), "bs", 0)
		# self.ub.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_SPODM_0"]), "bs", 0)
		# # self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_SPODM_1"]), "bs", 0)
		# # self.ub.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_SPODM_1"]), "bs", 0)
		# self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_SPODMpp"]), "bs", 0)
		# self.ub.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_SPODMpp"]), "bs", 0)
		# self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_FHL_0_0"]), "bs", 0)
		# self.ub.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_FHL_0_0"]), "bs", 0)
		# # self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_FHL_0_1"]), "bs", 0)
		# # self.ub.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_FHL_0_1"]), "bs", 0)
		# self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_FHL_1_0"]), "bs", 0)
		# self.ub.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_FHL_1_0"]), "bs", 0)
		# # self.lb.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_FHL_1_1"]), "bs", 0)
		# # self.ub.valuesIs(self.flexTFbaModel.rxnIdxs(["rxn_FEIST_FHL_1_1"]), "bs", 0)

		# self.flexTFbaModel.v_lowerIs(idxs = self.flexTFbaModel.rxnGroup("lowerMutable").idxs(), values = self.lb.mergedValues(self.flexTFbaModel.rxnGroup("lowerMutable").idxs()))
		# self.flexTFbaModel.v_upperIs(idxs = self.flexTFbaModel.rxnGroup("upperMutable").idxs(), values = self.ub.mergedValues(self.flexTFbaModel.rxnGroup("upperMutable").idxs()))

		# self.flexTFbaModel.solution()[self.flexTFbaModel.rxnIdxs(["mediaEx_FEIST_EX_glc(e)"])]
		# self.flexTFbaModel.solution()[self.flexTFbaModel.rxnIdxs(["mediaEx_FEIST_EX_o2(e)"])]
		# self.flexTFbaModel.solution()[self.flexTFbaModel.rxnIdxs(["g_bio"])]

		# import ipdb; ipdb.set_trace()

		# self.fba = setupFeist()


	def calculateRequest(self):
		self.ppi.requestAll()
		self.nmps.requestAll()
		self.atpRecyclingReactants.requestAll()

		self.biomassMetabolites.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# Store NMP/PPI counts for recycling

		ppiCount = self.ppi.count()
		nmpCounts = self.nmps.counts()
		hCount = self.h.count()
		piCount = self.pi.count()
		h2oCount = self.h2o.count()

		self.ppi.countIs(0)
		self.nmps.countsIs(0)
		self.h.countIs(0)
		self.pi.countIs(0)
		self.h2o.countIs(0)

		# Solve for metabolic fluxes

		effectiveBiomassObjective = self._computeEffectiveBiomass()

		if ASSUME_OPTIMAL_GROWTH:
			deltaMetabolitesNew = effectiveBiomassObjective * (
				1e-3 * self.nAvogadro * self.initialDryMass
				* np.exp(np.log(2)/self.cellCycleLen * self.time())
				* (np.exp(np.log(2)/self.cellCycleLen * self.time()) - 1)
				)

		else:
			raise NotImplementedError("Need to implement FBA")
		
		self.biomassMetabolites.countsInc(deltaMetabolitesNew.astype(np.int64))

		# NTP recycling
		
		if ppiCount >= nmpCounts.sum():
			self.ntps.countsInc(nmpCounts)
			self.ppi.countIs(ppiCount - nmpCounts.sum())

		else:
			recycledNmps = np.int64(nmpCounts/nmpCounts.sum() * ppiCount)

			self.ntps.countsInc(recycledNmps)
			self.nmps.countsIs(nmpCounts - recycledNmps)
			self.ppi.countIs(ppiCount - recycledNmps.sum())

		self.h.countIs(hCount)
		self.h2o.countIs(h2oCount)
		self.pi.countIs(piCount)
		adpsToRecycle = np.min(
			self.atpRecyclingReactants.counts()
			)
		self.atpRecyclingReactants.countsDec(adpsToRecycle)
		self.atpRecyclingProducts.countsInc(adpsToRecycle)

		# Write out effective biomass
		self.writeToListener("EffectiveBiomassObjective", "effectiveBiomassObjective", effectiveBiomassObjective)


	def _computeEffectiveBiomass(self):
		# Compute effective biomass objective

		metaboliteMasses = self.biomassMws.sum(axis = 1) / self.nAvogadro

		biomassObjective = self.wildtypeBiomassReactionSS

		# Solve for AA production
		aaIdxs = self.aaIdxsInWildTypeBiomass
		deltaAaMass = (
			np.dot(
				metaboliteMasses[aaIdxs],
				biomassObjective[aaIdxs] * 1e-3 * self.nAvogadro * self.initialDryMass
				) * np.exp(np.log(2) / self.cellCycleLen * self.time()) *
			(np.exp(np.log(2) / self.cellCycleLen * self.timeStepSec) - 1.0)
			)

		initialMetaboliteCounts = self.biomassMetabolites.counts()

		nominalBiomassProductionCoeff = (
			deltaAaMass + np.dot(metaboliteMasses[aaIdxs], initialMetaboliteCounts[aaIdxs])
			) / np.dot(metaboliteMasses[aaIdxs], biomassObjective[aaIdxs])

		deltaAaNew = np.int64(
			nominalBiomassProductionCoeff * biomassObjective[aaIdxs] -
			initialMetaboliteCounts[aaIdxs]
			)


		deltaAaMass = np.dot(metaboliteMasses[aaIdxs], deltaAaNew)
		deltaAaMassNonzero = np.dot(metaboliteMasses[aaIdxs], np.fmax(deltaAaNew, 0))

		deltaAaNew = np.int64(
			deltaAaMass / deltaAaMassNonzero * np.fmax(deltaAaNew, 0)
			)

		# Solve for DNTP production
		dntpIdxs = self.dntpIdxsInWildTypeBiomass
		deltaDntpMass = (
			np.dot(
				metaboliteMasses[dntpIdxs],
				biomassObjective[dntpIdxs] * 1e-3 * self.nAvogadro * self.initialDryMass
				) * np.exp(np.log(2) / self.cellCycleLen * self.time()) *
			(np.exp(np.log(2) / self.cellCycleLen * self.timeStepSec) - 1.0)
			)

		initialMetaboliteCounts = self.biomassMetabolites.counts()

		nominalBiomassProductionCoeff = (
			deltaDntpMass + np.dot(metaboliteMasses[dntpIdxs], initialMetaboliteCounts[dntpIdxs])
			) / np.dot(metaboliteMasses[dntpIdxs], biomassObjective[dntpIdxs])

		deltaDntpsNew = np.int64(
			nominalBiomassProductionCoeff * biomassObjective[dntpIdxs] -
			initialMetaboliteCounts[dntpIdxs]
			)


		deltaDntpsMass = np.dot(metaboliteMasses[dntpIdxs], deltaDntpsNew)
		deltaDntpsMassNonzero = np.dot(metaboliteMasses[dntpIdxs], np.fmax(deltaDntpsNew, 0))

		deltaDntpsNew = np.int64(
			deltaDntpsMass / deltaDntpsMassNonzero * np.fmax(deltaDntpsNew, 0)
			)

		# Solve for NTP production
		ntpIdxs = self.ntpIdxsInWildTypeBiomass
		deltaNtpMass = (
			np.dot(
				metaboliteMasses[ntpIdxs],
				biomassObjective[ntpIdxs] * 1e-3 * self.nAvogadro * self.initialDryMass
				) * np.exp(np.log(2) / self.cellCycleLen * self.time()) *
			(np.exp(np.log(2) / self.cellCycleLen * self.timeStepSec) - 1.0)
			)

		initialMetaboliteCounts = self.biomassMetabolites.counts()

		nominalBiomassProductionCoeff = (
			deltaNtpMass + np.dot(metaboliteMasses[ntpIdxs], initialMetaboliteCounts[ntpIdxs])
			) / np.dot(metaboliteMasses[ntpIdxs], biomassObjective[ntpIdxs])

		deltaNtpsNew = np.int64(
			nominalBiomassProductionCoeff * biomassObjective[ntpIdxs] -
			initialMetaboliteCounts[ntpIdxs]
			)


		deltaNtpsMass = np.dot(metaboliteMasses[ntpIdxs], deltaNtpsNew)
		deltaNtpsMassNonzero = np.dot(metaboliteMasses[ntpIdxs], np.fmax(deltaNtpsNew, 0))

		deltaNtpsNew = np.int64(
			deltaNtpsMass / deltaNtpsMassNonzero * np.fmax(deltaNtpsNew, 0)
			)

		# Solve for other metabolite production
		otherIdxs = self.notAADntpNtpIdxsInWildTypeBiomass
		deltaOtherMass = (
			np.dot(
				metaboliteMasses[otherIdxs],
				biomassObjective[otherIdxs] * 1e-3 * self.nAvogadro * self.initialDryMass
				) * np.exp(np.log(2) / self.cellCycleLen * self.time()) *
			(np.exp(np.log(2) / self.cellCycleLen * self.timeStepSec) - 1.0)
			)

		initialMetaboliteCounts = self.biomassMetabolites.counts()

		nominalBiomassProductionCoeff = (
			deltaOtherMass + np.dot(metaboliteMasses[otherIdxs], initialMetaboliteCounts[otherIdxs])
			) / np.dot(metaboliteMasses[otherIdxs], biomassObjective[otherIdxs])

		deltaOtherNew = np.int64(
			nominalBiomassProductionCoeff * biomassObjective[otherIdxs] -
			initialMetaboliteCounts[otherIdxs])


		deltaOtherMass = np.dot(metaboliteMasses[otherIdxs], deltaOtherNew)
		deltaOtherMassNonzero = np.dot(metaboliteMasses[otherIdxs], np.fmax(deltaOtherNew, 0))

		deltaOtherNew = np.int64(
			deltaOtherMass / deltaOtherMassNonzero * np.fmax(deltaOtherNew, 0)
			)

		# # Solve for non-DNTP metabolite production
		# notDntpIdxs = self.notDntpIdxsInWildTypeBiomass
		# deltaNotDntpMass = (
		# 	np.dot(
		# 		metaboliteMasses[notDntpIdxs],
		# 		biomassObjective[notDntpIdxs] * 1e-3 * self.nAvogadro * self.initialDryMass
		# 		) * np.exp(np.log(2) / self.cellCycleLen * self.time()) *
		# 	(np.exp(np.log(2) / self.cellCycleLen * self.timeStepSec) - 1.0)
		# 	)

		# initialMetaboliteCounts = self.biomassMetabolites.counts()

		# nominalBiomassProductionCoeff = (
		# 	deltaNotDntpMass + np.dot(metaboliteMasses[notDntpIdxs], initialMetaboliteCounts[notDntpIdxs])
		# 	) / np.dot(metaboliteMasses[notDntpIdxs], biomassObjective[notDntpIdxs])

		# deltaNotDntpsNew = np.int64(
		# 	nominalBiomassProductionCoeff * biomassObjective[notDntpIdxs] -
		# 	initialMetaboliteCounts[notDntpIdxs])


		# deltaNotDntpsMass = np.dot(metaboliteMasses[notDntpIdxs], deltaNotDntpsNew)
		# deltaNotDntpsMassNonzero = np.dot(metaboliteMasses[notDntpIdxs], np.fmax(deltaNotDntpsNew, 0))

		# deltaNotDntpsNew = np.int64(
		# 	deltaNotDntpsMass / deltaNotDntpsMassNonzero * np.fmax(deltaNotDntpsNew, 0)
		# 	)

		# Pools
		deltaPoolMetabolites = (
			(self.wildtypeBiomassPoolIncreases * 1e-3 * self.nAvogadro * self.initialDryMass) *
			np.exp(np.log(2) / self.cellCycleLen * self.time()) *
			(np.exp(np.log(2) / self.cellCycleLen * self.timeStepSec) - 1.0)
			)

		# Expected production
		expectedDeltaMetabolitesNew = (
			(self.wildtypeBiomassReactionSS * 1e-3 * self.nAvogadro * self.initialDryMass) *
			np.exp(np.log(2) / self.cellCycleLen * self.time()) *
			(np.exp(np.log(2) / self.cellCycleLen * self.timeStepSec) - 1.0)
			)

		# Aggregate all metabolite production (DNTPs and non-DNTPs)
		deltaMetabolitesNew = np.zeros_like(biomassObjective)
		deltaMetabolitesNew[aaIdxs] = deltaAaNew
		deltaMetabolitesNew[dntpIdxs] = deltaDntpsNew
		deltaMetabolitesNew[ntpIdxs] = deltaNtpsNew
		deltaMetabolitesNew[otherIdxs] = deltaOtherNew
		# deltaMetabolitesNew[notDntpIdxs] = deltaNotDntpsNew

		deltaMetabolitesNew += deltaPoolMetabolites

		effectiveBiomassObjective = deltaMetabolitesNew / (
			1e-3 * self.nAvogadro * self.initialDryMass
			* np.exp(np.log(2)/self.cellCycleLen * self.time())
			* (np.exp(np.log(2)/self.cellCycleLen * self.time()) - 1)
			)

		return effectiveBiomassObjective

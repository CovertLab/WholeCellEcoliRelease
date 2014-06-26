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
from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM

GAIN = 10

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

		effectiveBiomassObjective = deltaMetabolitesNew / (
			1e-3 * self.nAvogadro * self.initialDryMass
			* np.exp(np.log(2)/self.cellCycleLen * self.time())
			* (np.exp(np.log(2)/self.cellCycleLen * self.time()) - 1)
			)

		self.writeToListener("EffectiveBiomassObjective", "effectiveBiomassObjective", effectiveBiomassObjective)


	def _computeBiomassFromRequests(self):
		"""
		Computes the biomass objective dynamically as a function of the counts
		requested by processes, the counts remaining after the last time step,
		and the steady-state biomass objective.
		"""

		biomassObjective = self.wildtypeBiomassReactionSS.copy()

		requests = self.readFromListener("MetabolicDemands", "metaboliteRequests").sum(1)

		# print "Before: %0.10f" % (np.dot(self.biomassMws.sum(axis = 1) / 1000, self.wildtypeBiomassReaction))

		### Expected production

		## AAs
		expectedAAs = (
			np.round(
				self.wildtypeBiomassReactionSS[self.aaIdxsInWildTypeBiomass] * 1e-3 *
				self.nAvogadro * self.initialDryMass
				) *
			np.exp(np.log(2) / self.cellCycleLen * self.time()) *
			(np.exp(np.log(2) / self.cellCycleLen) - 1.0)
			)

		aasExpectedMass = np.dot(expectedAAs / self.nAvogadro, self.aaMws)

		## NTPs
		expectedNTPs = (
			np.round(
				self.wildtypeBiomassReactionSS[self.ntpIdxsInWildTypeBiomass] * 1e-3 *
				self.nAvogadro * self.initialDryMass
				) *
			np.exp(np.log(2) / self.cellCycleLen * self.time()) *
			(np.exp(np.log(2) / self.cellCycleLen) - 1.0)
			)

		ntpsExpectedMass = np.dot(expectedNTPs / self.nAvogadro, self.ntpMws)


		## dNTPs
		expectedDNTPs = (
			np.round(
				self.wildtypeBiomassReactionSS[self.dntpIdxsInWildTypeBiomass] * 1e-3 *
				self.nAvogadro * self.initialDryMass
				) *
			np.exp(np.log(2) / self.cellCycleLen * self.time()) *
			(np.exp(np.log(2) / self.cellCycleLen) - 1.0)
			)

		dntpsExpectedMass = np.dot(expectedDNTPs / self.nAvogadro, self.dntpMws)
		dntpsRequestedMass = np.dot(requests[self.dntpIdxsInWildTypeBiomass] / self.nAvogadro, self.dntpMws)

		dryMassExpectedProduction = (
			np.dot(
				self.biomassMws.sum(axis = 1) / 1000,
				(
					self.wildtypeBiomassReactionSS *
					self.initialDryMass *
					np.exp(np.log(2) / self.cellCycleLen * self.time()) *
					(np.exp(np.log(2) / self.cellCycleLen) - 1.0)
				)
				)
			)

		# Allocate masses
		# TODO: Try messing with NTP allocation
		dntpsAllocatedMass = dntpsRequestedMass
		aasAllocatedMass = aasExpectedMass + dntpsExpectedMass - dntpsAllocatedMass

		### For AAs

		deficitAAs = np.fmax(
			requests[self.aaIdxsInWildTypeBiomass] - self.aas._totalCount,
			0
			)

		np.set_printoptions(suppress = True)
		relativeAArequests = normalize(
			GAIN * deficitAAs + expectedAAs
			)

		biomassObjective[self.aaIdxsInWildTypeBiomass] = (
			countsFromMassAndExpression(
				(
					aasAllocatedMass / dryMassExpectedProduction *
					np.dot(
						self.biomassMws.sum(axis = 1) / 1000,
						self.wildtypeBiomassReactionSS
						)
				),
				self.aaMws,
				relativeAArequests,
				self.nAvogadro
				) *
			relativeAArequests *
			1000 / self.nAvogadro
			)

		### For NTPs

		deficitNTPs = np.fmax(
			requests[self.ntpIdxsInWildTypeBiomass] - self.nmps.counts() - self.ntps._totalCount,
			0
			)

		relativeNTPrequests = normalize(
			GAIN * deficitNTPs + expectedNTPs
			)

		biomassObjective[self.ntpIdxsInWildTypeBiomass] = (
			countsFromMassAndExpression(
				np.dot(
					self.ntpMws / 1000,
					self.wildtypeBiomassReactionSS[self.ntpIdxsInWildTypeBiomass]
					),
				self.ntpMws,
				relativeNTPrequests,
				self.nAvogadro
				) *
			relativeNTPrequests *
			1000 / self.nAvogadro
			)

		### For dNTPs

		deficitDNTPs = np.fmax(
			requests[self.dntpIdxsInWildTypeBiomass] - self.dntps._totalCount,
			0
			)

		relativeDNTPrequests = normalize(
			GAIN * deficitDNTPs + expectedDNTPs
			)

		biomassObjective[self.dntpIdxsInWildTypeBiomass] = (
			countsFromMassAndExpression(
				(
					dntpsAllocatedMass / dryMassExpectedProduction *
					np.dot(
						self.biomassMws.sum(axis = 1) / 1000,
						self.wildtypeBiomassReactionSS
						)
				),
				self.dntpMws,
				relativeDNTPrequests,
				self.nAvogadro
				) *
			relativeDNTPrequests *
			1000 / self.nAvogadro
			)

		# print "After: %0.10f" % (np.dot(self.biomassMws.sum(axis = 1) / 1000, self.wildtypeBiomassReaction))

		return biomassObjective


	def _computeBiomassFromLeftovers(self):
		"""
		Computes the biomass objective dynamically based on the counts 
		remaining from the last time step and the steady-state biomass 
		objective.
		"""

		biomassObjective = self.wildtypeBiomassReactionSS.copy()

		# print "Before: %0.10f" % (np.dot(self.biomassMws.sum(axis = 1) / 1000, self.wildtypeBiomassReaction))


		### For AAs

		expectedAAs = (
			np.round(
				self.wildtypeBiomassReactionSS[self.aaIdxsInWildTypeBiomass] * 1e-3 *
				self.nAvogadro * self.initialDryMass
				) *
			np.exp(np.log(2) / self.cellCycleLen * self.time()) *
			(np.exp(np.log(2) / self.cellCycleLen) - 1.0)
			)

		deficitAAs = np.fmax(
			2*expectedAAs - self.aas._totalCount,
			0
			)

		np.set_printoptions(suppress = True)
		relativeAArequests = normalize(
			GAIN * deficitAAs + expectedAAs
			)

		biomassObjective[self.aaIdxsInWildTypeBiomass] = (
			countsFromMassAndExpression(
				np.dot(
					self.aaMws / 1000,
					self.wildtypeBiomassReactionSS[self.aaIdxsInWildTypeBiomass]
					),
				self.aaMws,
				relativeAArequests,
				self.nAvogadro
				) *
			relativeAArequests *
			1000 / self.nAvogadro
			)

		### For NTPs

		expectedNTPs = (
			np.round(
				self.wildtypeBiomassReactionSS[self.ntpIdxsInWildTypeBiomass] * 1e-3 *
				self.nAvogadro * self.initialDryMass
				) *
			np.exp(np.log(2) / self.cellCycleLen * self.time()) *
			(np.exp(np.log(2) / self.cellCycleLen) - 1.0)
			)

		deficitNTPs = np.fmax(
			2*expectedNTPs - self.nmps.counts() - self.ntps._totalCount,
			0
			)

		relativeNTPrequests = normalize(
			GAIN * deficitNTPs + expectedNTPs
			)

		biomassObjective[self.ntpIdxsInWildTypeBiomass] = (
			countsFromMassAndExpression(
				np.dot(
					self.ntpMws / 1000,
					self.wildtypeBiomassReactionSS[self.ntpIdxsInWildTypeBiomass]
					),
				self.ntpMws,
				relativeNTPrequests,
				self.nAvogadro
				) *
			relativeNTPrequests *
			1000 / self.nAvogadro
			)

		### For dNTPs

		expectedDNTPs = (
			np.round(
				self.wildtypeBiomassReactionSS[self.dntpIdxsInWildTypeBiomass] * 1e-3 *
				self.nAvogadro * self.initialDryMass
				) *
			np.exp(np.log(2) / self.cellCycleLen * self.time()) *
			(np.exp(np.log(2) / self.cellCycleLen) - 1.0)
			)

		deficitDNTPs = np.fmax(
			2*expectedDNTPs - self.dntps._totalCount,
			0
			)

		relativeDNTPrequests = normalize(
			GAIN * deficitDNTPs + expectedDNTPs
			)

		biomassObjective[self.dntpIdxsInWildTypeBiomass] = (
			countsFromMassAndExpression(
				np.dot(
					self.dntpMws / 1000,
					self.wildtypeBiomassReactionSS[self.dntpIdxsInWildTypeBiomass]
					),
				self.dntpMws,
				relativeDNTPrequests,
				self.nAvogadro
				) *
			relativeDNTPrequests *
			1000 / self.nAvogadro
			)

		# print "After: %0.10f" % (np.dot(self.biomassMws.sum(axis = 1) / 1000, self.wildtypeBiomassReaction))


		return biomassObjective

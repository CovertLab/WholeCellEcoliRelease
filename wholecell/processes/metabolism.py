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

from itertools import izip

import numpy as np

import wholecell.processes.process

from wholecell.utils.random import stochasticRound
from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM

from wholecell.utils.modular_fba import FluxBalanceAnalysis

class Metabolism(wholecell.processes.process.Process):
	""" Metabolism """

	_name = "Metabolism"

	# Constructor
	def __init__(self):

		super(Metabolism, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Metabolism, self).initialize(sim, kb)

		# Load constants
		self.nAvogadro = kb.nAvogadro.to('1 / mole').magnitude
		self.cellDensity = kb.cellDensity.to("g/L").magnitude
		
		self.metabolitePoolIDs = kb.metabolitePoolIDs
		self.targetConcentrations = kb.metabolitePoolConcentrations.to("mole/L").magnitude
		
		# Set up FBA solver

		import re
		rxns = [
			x for x in kb.metabolismBiochemicalReactions
			if (
				not re.match(".*_[0-9]$", x["id"])
				or x["id"].endswith("_0")
				or "PFK_2" in x["id"]
				)
			]

		mediaEx = kb.metabolismMediaEx

		reactionStoich = {
			rxn["id"]:
			{entry["molecule"]+"["+entry["location"]+"]":entry["coeff"] for entry in rxn["stoichiometry"]}
			for rxn in rxns
			if len(rxn["stoichiometry"]) > 1 # no exchange reactions!
			}

		reversibleReactions = [rxn["id"] for rxn in rxns if not rxn["dir"]]

		externalExchangedMolecules = [rxn["met"] for rxn in mediaEx]

		objective = {
			moleculeID:coeff for moleculeID, coeff in
			izip(self.metabolitePoolIDs, self.targetConcentrations)
			}

		reactionEnzymes = {
			reactionID:enzymeID
			for reactionID, enzymeID in izip(kb.metabolismReactionIds, kb.metabolismReactionEnzymes)
			if (enzymeID is not None) and reactionStoich.has_key(reactionID)
			}

		dt = self.timeStepSec # rates are per-second but media exchange fluxes are per-hour
		reactionRates = {
			reactionID:rate * dt
			for reactionID, rate in izip(kb.metabolismReactionIds, kb.metabolismReactionKcat)
			if reactionStoich.has_key(reactionID) and rate > 0
			}

		masses = kb.getMass(externalExchangedMolecules).to("gram/mole").magnitude

		moleculeMasses = {moleculeID:masses[index]
			for index, moleculeID in enumerate(externalExchangedMolecules)}

		self.fba = FluxBalanceAnalysis(
			reactionStoich,
			externalExchangedMolecules,
			objective,
			objectiveType = "pools",
			reversibleReactions = reversibleReactions,
			reactionEnzymes = reactionEnzymes,
			reactionRates = reactionRates,
			moleculeMasses = moleculeMasses
			)

		# Set constraints
		## External molecules
		externalMoleculeIDs = self.fba.externalMoleculeIDs()

		unconstrainedExchange = ("CA2[e]", "CL[e]", "CO2[e]", "COBALT2[e]", 
			"CU2[e]", "FE2[e]", "FE3[e]", "H[e]", "H2O[e]", "K[e]", "MG2[e]",
			"MN2[e]", "MOBD[e]", "NA1[e]", "NH4[e]", "PI[e]", "SO4[e]", "TUNGS[e]",
			"ZN2[e]")

		# (flux) mmol/gDCW/hr
		# * 1 hr / 3600 s = mmol/gDCW/s
		# * 1 mol / 1000 mmol = mol/gDCW/s
		# * dt = mol/gDCW
		# * initDry/initTotal = mol/g
		# * cellDensity = mol/L = M

		initWaterMass = kb.avgCellWaterMassInit.to('gram * water_gram / DCW_gram').magnitude
		initDryMass = kb.avgCellDryMassInit.to('gram').magnitude

		initCellMass = (
			initWaterMass
			+ initDryMass
			)

		coeff = 1/3600 * 1e-3 * dt * initDryMass / initCellMass * self.cellDensity

		self._coeff = coeff

		constrainedExchange = {
			"CBL1[e]":0.01 * coeff,
			"GLC-D[e]":8 * coeff,
			"O2[e]":18.5 * coeff
			}

		externalMoleculeLevels = np.zeros(len(externalMoleculeIDs), np.float64)

		for index, moleculeID in enumerate(externalMoleculeIDs):
			if moleculeID in unconstrainedExchange:
				externalMoleculeLevels[index] = np.inf

			elif constrainedExchange.has_key(moleculeID):
				externalMoleculeLevels[index] = constrainedExchange[moleculeID]

		self.fba.externalMoleculeLevelsIs(externalMoleculeLevels)

		## Set Feist's forced reactions

		### ATP maintenance
		self.fba.minReactionFluxIs("FEIST_ATPM", 8.39 * coeff)
		self.fba.maxReactionFluxIs("FEIST_ATPM", 8.39 * coeff)

		### Arbitrarily disabled reactions
		disabledReactions = ("FEIST_CAT_0", "FEIST_SPODM_0", "FEIST_SPODMpp", "FEIST_FHL_0_0", "FEIST_FHL_1_0")

		for reactionID in disabledReactions:
			self.fba.maxReactionFluxIs(reactionID, 0)

		## Set enzymes unlimited
		self.fba.enzymeLevelsIs(np.inf)

		# Views
		self.metabolites = self.bulkMoleculesView(self.fba.outputMoleculeIDs())
		self.poolMetabolites = self.bulkMoleculesView(self.metabolitePoolIDs)

		assert self.fba.outputMoleculeIDs() == self.fba.internalMoleculeIDs()

		# Set the priority to a low value

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)


	def calculateRequest(self):
		self.metabolites.requestAll()


	# Calculate temporal evolution
	def evolveState(self):
		# Solve for metabolic fluxes

		metaboliteCountsInit = self.metabolites.counts()
		poolCounts = self.poolMetabolites.counts()

		cellMass = self.readFromListener("Mass", "cellMass") * 1e-15 # fg to g

		cellVolume = cellMass / self.cellDensity

		# if self.time() < 2:
		# 	print poolCounts / self.nAvogadro / cellVolume

		countsToMolar = 1 / (self.nAvogadro * cellVolume)

		self.fba.internalMoleculeLevelsIs(
			metaboliteCountsInit * countsToMolar
			)
			
		self.fba.run()

		deltaMetabolites = self.fba.outputMoleculeLevelsChange() / countsToMolar

		metaboliteCountsFinal = np.fmax(stochasticRound(
			self.randomState,
			metaboliteCountsInit + deltaMetabolites
			), 0).astype(np.int64)

		self.metabolites.countsIs(metaboliteCountsFinal)

		# print "mass: {:0.2f}".format(self.fba.massAccumulated()/self._coeff)
		# print "glucose: {:0.2f}".format(self.fba.externalExchangeFlux("GLC-D[e]")/self._coeff)
		# print "oxygen: {:0.2f}".format(self.fba.externalExchangeFlux("O2[e]")/self._coeff)
		# print "cbl1: {:0.2f}".format(self.fba.externalExchangeFlux("CBL1[e]")/self._coeff)

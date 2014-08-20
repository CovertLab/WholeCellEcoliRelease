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
from wholecell.utils import units

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
		self.nAvogadro = kb.nAvogadro.asNumber(1 / units.mol)
		self.cellDensity = kb.cellDensity.asNumber(units.g/units.L)
		
		self.metabolitePoolIDs = kb.metabolitePoolIDs
		self.targetConcentrations = kb.metabolitePoolConcentrations.asNumber(units.mol/units.L)
		
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

		masses = kb.getMass(externalExchangedMolecules).asNumber(units.g/units.mol)

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
			"ZN2[e]", "SELNP[e]",) #"INOST[e]", "23CAMP[e]", )

		# (flux) mmol/gDCW/hr
		# * 1 hr / 3600 s = mmol/gDCW/s
		# * 1 mol / 1000 mmol = mol/gDCW/s
		# * dt = mol/gDCW
		# * initDry/initTotal = mol/g
		# * cellDensity = mol/L = M

		initWaterMass = kb.avgCellWaterMassInit.asNumber(units.g)
		initDryMass = kb.avgCellDryMassInit.asNumber(units.g)

		initCellMass = (
			initWaterMass
			+ initDryMass
			)

		coeff = 1/3600 * 1e-3 * dt * initDryMass / initCellMass * self.cellDensity

		# self._coeff = coeff

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
		# NOTE: all maintenance is handled in the AtpUsage process
		self.fba.maxReactionFluxIs("FEIST_ATPM", 0)

		### Arbitrarily disabled reactions
		disabledReactions = ("FEIST_CAT_0", "FEIST_SPODM_0", "FEIST_SPODMpp",
			"FEIST_FHL_0_0", "FEIST_FHL_1_0")

		for reactionID in disabledReactions:
			self.fba.maxReactionFluxIs(reactionID, 0)

		## Set enzymes unlimited
		self.fba.enzymeLevelsIs(np.inf)

		# Views
		self.metabolites = self.bulkMoleculesView(self.fba.outputMoleculeIDs())
		self.poolMetabolites = self.bulkMoleculesView(self.metabolitePoolIDs)

		outputMoleculeIDs = self.fba.outputMoleculeIDs()

		assert outputMoleculeIDs == self.fba.internalMoleculeIDs()

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

		# TODO: report as reactions (#) per second & store volume elsewhere
		self.writeToListener("FBAResults", "reactionFluxes",
			self.fba.reactionFluxes() / self.timeStepSec)
		self.writeToListener("FBAResults", "externalExchangeFluxes",
			self.fba.externalExchangeFluxes() / self.timeStepSec)
		self.writeToListener("FBAResults", "objectiveValue",
			self.fba.objectiveValue() / deltaMetabolites.size) # divide to normalize by number of metabolites
		self.writeToListener("FBAResults", "outputFluxes",
			self.fba.outputMoleculeLevelsChange() / self.timeStepSec)

		# NOTE: the calculation for the objective components doesn't yet have
		# an interface, since it will vary in calculation and shape for every
		# objective type

		objectiveComponents_raw = (np.array(self.fba._f).flatten() * self.fba._solutionFluxes)[self.fba._objIndexes]
		objectiveComponents = objectiveComponents_raw[::2] + objectiveComponents_raw[1::2]

		self.writeToListener("FBAResults", "objectiveComponents",
			objectiveComponents
			)

		# TODO:
		# - which media exchanges/reactions are limiting, if any
		# - objective details (value, component values)

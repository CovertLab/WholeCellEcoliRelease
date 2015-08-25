#!/usr/bin/env python

"""
Metabolism

Metabolism sub-model. Encodes molecular simulation of microbial metabolism using flux-balance analysis.

TODO:
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
from wholecell.utils.enzymeKinetics import EnzymeKinetics

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g
USE_RATELIMITS = False # Enable/disable kinetic rate limits in the model

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
		self.nAvogadro = kb.constants.nAvogadro.asNumber(1 / COUNTS_UNITS)
		self.cellDensity = kb.constants.cellDensity.asNumber(MASS_UNITS/VOLUME_UNITS)

		self.metabolitePoolIDs = kb.process.metabolism.metabolitePoolIDs
		self.targetConcentrations = kb.process.metabolism.metabolitePoolConcentrations.asNumber(COUNTS_UNITS/VOLUME_UNITS)

		# Load enzyme kinetic rate information
		self.reactionRateInfo = kb.process.metabolism.reactionRateInfo
		self.enzymesWithKineticInfo = kb.process.metabolism.enzymesWithKineticInfo["enzymes"]
		self.constraintIDs = kb.process.metabolism.constraintIDs
		self.constraintToReactionDict = kb.process.metabolism.constraintToReactionDict

		objective = dict(zip(
			self.metabolitePoolIDs,
			self.targetConcentrations
			))

		# TODO: make kb method?
		extIDs = kb.process.metabolism.externalExchangeMolecules

		moleculeMasses = dict(zip(
			extIDs,
			kb.getter.getMass(extIDs).asNumber(MASS_UNITS/COUNTS_UNITS)
			))

		initWaterMass = kb.mass.avgCellWaterMassInit
		initDryMass = kb.mass.avgCellDryMassInit

		initCellMass = (
			initWaterMass
			+ initDryMass
			)

		energyCostPerWetMass = kb.constants.darkATP * initDryMass / initCellMass

		# Set up FBA solver
		self.fba = FluxBalanceAnalysis(
			kb.process.metabolism.reactionStoich.copy(), # TODO: copy in class
			kb.process.metabolism.externalExchangeMolecules,
			objective,
			objectiveType = "pools",
			reversibleReactions = kb.process.metabolism.reversibleReactions,
			moleculeMasses = moleculeMasses,
			# maintenanceCost = energyCostPerWetMass.asNumber(COUNTS_UNITS/MASS_UNITS), # mmol/gDCW TODO: get real number
			# maintenanceReaction = {
			# 	"ATP[c]":-1, "WATER[c]":-1, "ADP[c]":+1, "Pi[c]":+1
			# 	} # TODO: move to KB TODO: check reaction stoich
			)

		# Set up enzyme kinetics object
		self.enzymeKinetics = EnzymeKinetics(kb, reactionIDs = self.fba.reactionIDs(), metaboliteIDs = self.fba.outputMoleculeIDs(), kcatOnly=False)

		# Determine which kinetic limits to use
		self.reactionsWithKineticLimits = [True]*len(self.fba.reactionIDs())

		# halfPoint = int(np.floor(len(self.fba.reactionIDs())/2))
		# quarterPoint = int(np.floor(len(self.fba.reactionIDs())/4))
		# eigthPoint = int(np.floor(len(self.fba.reactionIDs())/8))
		# sixteenthPoint = int(np.floor(len(self.fba.reactionIDs())/16))
		# thirtysecondPoint = int(np.floor(len(self.fba.reactionIDs())/32))
		# sixtyfourthPoint = int(np.floor(len(self.fba.reactionIDs())/64))
		# a128thPoint = int(np.floor(len(self.fba.reactionIDs())/128))
		# a256thPoint = int(np.floor(len(self.fba.reactionIDs())/256))

		# start = 0 
		# end = start + 3*eigthPoint

		# for index in xrange(start,end):
		# 	self.reactionsWithKineticLimits[index] = True

		# Set constraints
		## External molecules
		externalMoleculeIDs = self.fba.externalMoleculeIDs()

		coefficient = initDryMass / initCellMass * kb.constants.cellDensity * (self.timeStepSec * units.s)

		externalMoleculeLevels = kb.process.metabolism.exchangeConstraints(
			externalMoleculeIDs,
			coefficient,
			COUNTS_UNITS / VOLUME_UNITS
			)

		self.fba.externalMoleculeLevelsIs(externalMoleculeLevels)

		## Set enzymes unlimited
		self.fba.enzymeLevelsIs(np.inf)

		# Views
		self.metaboliteNames = self.fba.outputMoleculeIDs()
		self.metabolites = self.bulkMoleculesView(self.metaboliteNames)
		self.poolMetabolites = self.bulkMoleculesView(self.metabolitePoolIDs)
		self.enzymeNames = self.enzymesWithKineticInfo
		self.enzymes = self.bulkMoleculesView(self.enzymeNames)

		outputMoleculeIDs = self.fba.outputMoleculeIDs()

		assert outputMoleculeIDs == self.fba.internalMoleculeIDs()

		# Set the priority to a low value

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)


	def calculateRequest(self):
		self.metabolites.requestAll()
		self.enzymes.requestAll()

	# Calculate temporal evolution
	def evolveState(self):

		# Solve for metabolic fluxes

		metaboliteCountsInit = self.metabolites.counts()

		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg).asNumber(MASS_UNITS)

		cellVolume = cellMass / self.cellDensity

		countsToMolar = 1 / (self.nAvogadro * cellVolume)


		#  Find metabolite concentrations from metabolite counts

		metaboliteConcentrations = metaboliteCountsInit * countsToMolar

		self.fba.internalMoleculeLevelsIs(
			metaboliteConcentrations
			)

		#  Find enzyme concentrations from enzyme counts

		enzymeCountsInit = self.enzymes.counts()

		enzymeConcentrations = enzymeCountsInit * countsToMolar

		defaultRate = self.enzymeKinetics.defaultRate

		# Combine the enzyme concentrations, substrate concentrations, and the default rate into one vector
		inputConcentrations = np.concatenate((enzymeConcentrations,metaboliteConcentrations,[defaultRate]), axis=1)

		# Find reaction rate limits
		self.reactionRates = self.enzymeKinetics.rateFunction(*inputConcentrations)

		# Find rate limits for all constraints
		self.allConstraintsLimits = self.enzymeKinetics.allRatesFunction(*inputConcentrations)

		# Find per-enzyme reaction rates
		perEnzymeInputConcentrations = np.concatenate((([1]*len(enzymeConcentrations)),metaboliteConcentrations,[defaultRate]), axis=1)
		self.perEnzymeRates = self.enzymeKinetics.rateFunction(*perEnzymeInputConcentrations)

		# Set the rate limits only if the option flag is enabled
		if USE_RATELIMITS:
			# Set max reaction fluxes for enzymes for which kinetics are known
			for index, reactionID in enumerate(self.fba.reactionIDs()):
				# Only use this kinetic limit if it's enabled
				if self.reactionsWithKineticLimits[index]:
					# Make sure to never set negative maximum rates
					assert (self.reactionRates[0][index] >= 0 and self.reactionRates[0][index] != np.nan)
					# Set the max reaction rate for this reaction
					self.fba.maxReactionFluxIs(reactionID, self.reactionRates[0][index], raiseForReversible = False)
				else:
					self.fba.maxReactionFluxIs(reactionID, defaultRate, raiseForReversible = False)


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
		# self.writeToListener("FBAResults", "objectiveValue", # TODO
		# 	self.fba.objectiveValue() / deltaMetabolites.size) # divide to normalize by number of metabolites
		self.writeToListener("FBAResults", "outputFluxes",
			self.fba.outputMoleculeLevelsChange() / self.timeStepSec)

		self.writeToListener("FBAResults", "outputFluxes",
			self.fba.outputMoleculeLevelsChange() / self.timeStepSec)

		self.writeToListener("EnzymeKinetics", "reactionRates",
			self.reactionRates)

		self.writeToListener("EnzymeKinetics", "perEnzymeRates",
			self.perEnzymeRates)

		self.writeToListener("EnzymeKinetics", "allConstraintsLimits",
			self.allConstraintsLimits)

		self.writeToListener("EnzymeKinetics", "metaboliteCountsInit",
			metaboliteCountsInit)

		self.writeToListener("EnzymeKinetics", "metaboliteConcentrations",
			metaboliteConcentrations)

		self.writeToListener("EnzymeKinetics", "countsToMolar",
			countsToMolar)




		# TODO
		# NOTE: the calculation for the objective components doesn't yet have
		# an interface, since it will vary in calculation and shape for every
		# objective type

		# objectiveComponents_raw = (np.array(self.fba._f).flatten() * self.fba._solutionFluxes)[self.fba._objIndexes]
		# objectiveComponents = objectiveComponents_raw[::2] + objectiveComponents_raw[1::2]

		# self.writeToListener("FBAResults", "objectiveComponents",
		# 	objectiveComponents
		# 	)

		# TODO:
		# - which media exchanges/reactions are limiting, if any
		# - objective details (value, component values)

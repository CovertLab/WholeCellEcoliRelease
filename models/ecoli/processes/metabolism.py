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

from wholecell.utils.enzymeKinetics import *

COUNTS_UNITS = units.mmol
VOLUME_UNITS = units.L
MASS_UNITS = units.g

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



		objective = {
			moleculeID:coeff for moleculeID, coeff in
			izip(self.metabolitePoolIDs, self.targetConcentrations)
			}

		# Set up FBA solver

		self.fba = FluxBalanceAnalysis(
			kb.process.metabolism.reactionStoich.copy(), # TODO: copy in class
			kb.process.metabolism.externalExchangeMolecules,
			objective,
			objectiveType = "pools",
			reversibleReactions = kb.process.metabolism.reversibleReactions,
			# reactionEnzymes = kb.process.metabolism.reactionEnzymes.copy(), # TODO: copy in class
			# reactionRates = kb.process.metabolism.reactionRates(self.timeStepSec * units.s),
			# moleculeMasses = kb.process.metabolism.exchangeMasses(MASS_UNITS / COUNTS_UNITS)
			)

		# Set constraints
		## External molecules
		externalMoleculeIDs = self.fba.externalMoleculeIDs()

		initWaterMass = kb.mass.avgCellWaterMassInit
		initDryMass = kb.mass.avgCellDryMassInit

		initCellMass = (
			initWaterMass
			+ initDryMass
			)

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
		self.metabolites = self.bulkMoleculesView(self.fba.outputMoleculeIDs())
		self.metaboliteNames = self.fba.outputMoleculeIDs()
		self.poolMetabolites = self.bulkMoleculesView(self.metabolitePoolIDs)
		self.enzymes = self.bulkMoleculesView(self.enzymesWithKineticInfo)

		outputMoleculeIDs = self.fba.outputMoleculeIDs()

		assert outputMoleculeIDs == self.fba.internalMoleculeIDs()

		# Set the priority to a low value

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)

		# Make a dictionary mapping a substrate ID to it's index in self.metabolites()
		self.metaboliteIndexDict = {}
		for index,name in enumerate(self.metaboliteNames):
			self.metaboliteIndexDict[name] = index

		# Make a dictionary mapping an enzyme ID to it's index in self.enzymes()
		self.enzymeIndexDict = {}
		for index,name in enumerate(self.enzymesWithKineticInfo):
			self.enzymeIndexDict[name] = index

	def calculateRequest(self):
		self.metabolites.requestAll()
		self.enzymes.requestAll()

	# Calculate temporal evolution
	def evolveState(self):

		# Solve for metabolic fluxes

		metaboliteCountsInit = self.metabolites.counts()
		poolCounts = self.poolMetabolites.counts()

		cellMass = (self.readFromListener("Mass", "cellMass") * units.fg).asNumber(MASS_UNITS)

		cellVolume = cellMass / self.cellDensity

		# if self.time() < 2:
		# 	print poolCounts / self.nAvogadro / cellVolume

		countsToMolar = 1 / (self.nAvogadro * cellVolume)


		#  Find metabolite concentrations from metabolite counts

		metaboliteConcentrations = metaboliteCountsInit * countsToMolar

		self.fba.internalMoleculeLevelsIs(
			metaboliteConcentrations
			)

		#  Find enzyme concentrations from enzyme counts

		enzymeCountsInit = self.enzymes.counts()

		enzymeConcentrations = enzymeCountsInit * countsToMolar


		# find reaction rate limits

		# This list will hold the limits - default limit to infinity
		self.reactionRates = np.ones(len(self.fba.reactionIDs()))*np.inf

		for index, reactionID in enumerate(self.fba.reactionIDs()):
			rateInfo = {}
			try:
				rateInfo = self.reactionRateInfo[reactionID]
			except:
				continue

			# enzymes, substrates = enzymeRateValuesNeededLookup(rateInfo)

			substrateIDs = rateInfo["substrateIDs"]
			substrateIXs = [self.metaboliteIndexDict[x] for x in substrateIDs]
			substrateConcArray = [metaboliteConcentrations[i] for i in substrateIXs]

			enzymeIDs = rateInfo["enzymeIDs"]
			enzymeIXs = [self.enzymeIndexDict[x] for x in enzymeIDs]
			enzymeConcArray = [enzymeConcentrations[i] for i in enzymeIXs]

			rate = enzymeRate(rateInfo, enzymeConcArray, substrateConcArray)

			self.reactionRates[index] = rate




		# self.fba.run()

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
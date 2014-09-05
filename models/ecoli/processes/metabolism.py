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

		objective = {
			moleculeID:coeff for moleculeID, coeff in
			izip(self.metabolitePoolIDs, self.targetConcentrations)
			}
		
		# Set up FBA solver

		self.fba = FluxBalanceAnalysis(
			kb.metabolismReactionStoich.copy(), # TODO: copy in class
			kb.metabolismExternalExchangeMolecules,
			objective,
			objectiveType = "pools",
			reversibleReactions = kb.metabolismReversibleReactions,
			# reactionEnzymes = kb.metabolismReactionEnzymes.copy(), # TODO: copy in class
			# reactionRates = kb.metabolismReactionRates(self.timeStepSec * units.s),
			# moleculeMasses = kb.metabolismExchangeMasses(units.g / units.mol)
			)

		# Set constraints
		## External molecules
		externalMoleculeIDs = self.fba.externalMoleculeIDs()

		initWaterMass = kb.avgCellWaterMassInit
		initDryMass = kb.avgCellDryMassInit

		initCellMass = (
			initWaterMass
			+ initDryMass
			)

		coefficient = initDryMass / initCellMass * kb.cellDensity * (self.timeStepSec * units.s)

		externalMoleculeLevels = kb.metabolismExchangeConstraints(
			externalMoleculeIDs,
			coefficient,
			units.mol / units.L
			)

		self.fba.externalMoleculeLevelsIs(externalMoleculeLevels)

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

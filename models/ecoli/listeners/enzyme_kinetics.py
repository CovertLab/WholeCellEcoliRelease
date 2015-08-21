#!/usr/bin/env python

"""
EnzymeKinetics

EnzymeKinetics listener. Tracks information about enzyme kinetics.

@author: Morgan Paull
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/02/2015
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener
from wholecell.utils.fitting import normalize
from wholecell.utils import units

class EnzymeKinetics(wholecell.listeners.listener.Listener):
	""" EnzymeKinetics """

	_name = 'EnzymeKinetics'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(EnzymeKinetics, self).__init__(*args, **kwargs)

	# Construct object graph
	def initialize(self, sim, kb):
		super(EnzymeKinetics, self).initialize(sim, kb)

		self.metabolism = sim.processes["Metabolism"]
		self.reactionRateInfo = kb.process.metabolism.reactionRateInfo
		self.metaboliteIDs = kb.process.metabolism.metabolitePoolIDs
		self.reactionStoich = kb.process.metabolism.reactionStoich.copy()
		self.externalExchangeMolecules = kb.process.metabolism.externalExchangeMolecules
		self.metabolitePoolIDs = kb.process.metabolism.metabolitePoolIDs
		self.targetConcentrations = kb.process.metabolism.metabolitePoolConcentrations.asNumber(units.mmol/units.L)
		self.constraintToReactionDict = kb.process.metabolism.constraintToReactionDict

		self.objective = dict(zip(
			kb.process.metabolism.metabolitePoolIDs,
			self.targetConcentrations
			))

		self.reversibleReactions = kb.process.metabolism.reversibleReactions

		extIDs = kb.process.metabolism.externalExchangeMolecules
		
		self.moleculeMasses = dict(zip(
			extIDs,
			kb.getter.getMass(extIDs).asNumber(units.g/units.mmol)
			))

	# Allocate memory
	# In case things are of unknown size, write them here
	# Dummy values for what will be writen to output table
	# prep variables with zeros or NaNs, with correct size
	# to be filled later
	def allocate(self):
		super(EnzymeKinetics, self).allocate()

		self.reactionRates = np.zeros(len(self.metabolism.fba.reactionIDs()), np.float64)
		self.allConstraintsLimits = np.zeros(len(self.reactionRateInfo), np.float64)
		self.perEnzymeRates = np.zeros(len(self.metabolism.fba.reactionIDs()), np.float64)
		self.reactionIDs = self.metabolism.fba.reactionIDs()
		self.constraintIDs = self.metabolism.constraintIDs
		self.metaboliteCountsInit = np.zeros(len(self.metaboliteIDs), np.float64)
		self.metaboliteConcentrations = np.zeros(len(self.metaboliteIDs), np.float64)
		self.countsToMolar = np.zeros(1, np.float64)


	def update(self):
		pass

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			reactionIDs = self.reactionIDs,
			constraintIDs = self.constraintIDs,
			reactionStoich = self.reactionStoich,
			externalExchangeMolecules = self.externalExchangeMolecules,
			objective = self.objective,
			reversibleReactions = self.reversibleReactions,
			moleculeMasses = self.moleculeMasses,
			metabolitePoolIDs = self.metabolitePoolIDs,
			targetConcentrations = self.targetConcentrations,
			constraintToReactionDict = self.constraintToReactionDict,
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			reactionRates = self.reactionRates,
			perEnzymeRates = self.perEnzymeRates,
			allConstraintsLimits = self.allConstraintsLimits,
			metaboliteCountsInit = self.metaboliteCountsInit,
			metaboliteConcentrations = self.metaboliteConcentrations,
			countsToMolar = self.countsToMolar,
			)
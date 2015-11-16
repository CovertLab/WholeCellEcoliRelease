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
		self.constraintToReactionDict = kb.process.metabolism.constraintToReactionDict

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
		self.metaboliteCountsFinal = np.zeros(len(self.metaboliteIDs), np.float64)
		self.metaboliteConcentrations = np.zeros(len(self.metaboliteIDs), np.float64)
		self.enzymeCountsInit = np.zeros(len(self.metabolism.enzymesWithKineticInfo), np.float64)

		self.countsToMolar = np.zeros(1, np.float64)
		self.counts_units = "                                          " # Placeholder string longer than any unit name
		self.mass_units = "                                          " # Placeholder string longer than any unit name
		self.volume_units = "                                          " # Placeholder string longer than any unit name


	def update(self):
		pass

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			reactionIDs = self.reactionIDs,
			constraintIDs = self.constraintIDs,
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
			metaboliteCountsFinal = self.metaboliteCountsFinal,
			metaboliteConcentrations = self.metaboliteConcentrations,
			countsToMolar = self.countsToMolar,
			counts_units = self.counts_units,
			mass_units = self.mass_units,
			volume_units = self.volume_units,
			enzymeCountsInit = self.enzymeCountsInit
			)
#!/usr/bin/env python

"""
FBAResults

Records dynamics of FBA output.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/24/2014
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

class FBAResults(wholecell.listeners.listener.Listener):
	""" FBAResults """

	_name = "FBAResults"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(FBAResults, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(FBAResults, self).initialize(sim, kb)

		self.metabolism = sim.processes["Metabolism"]

		self.objectiveValue = 0.0


	# Allocate memory
	def allocate(self):
		super(FBAResults, self).allocate()

		fba = self.metabolism.fba

		self.reactionIDs = fba.reactionIDs()
		self.externalMoleculeIDs = fba.externalMoleculeIDs()
		self.outputMoleculeIDs = fba.outputMoleculeIDs()

		self.reactionFluxes = np.zeros(len(self.reactionIDs), np.float64)
		self.externalExchangeFluxes = np.zeros(len(self.externalMoleculeIDs), np.float64)
		self.outputFluxes = np.zeros(len(self.outputMoleculeIDs), np.float64)
		self.objectiveComponents = np.zeros_like(self.outputFluxes)


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			reactionIDs = self.reactionIDs,
			externalMoleculeIDs = self.externalMoleculeIDs,
			outputMoleculeIDs = self.outputMoleculeIDs
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			reactionFluxes = self.reactionFluxes,
			externalExchangeFluxes = self.externalExchangeFluxes,
			outputFluxes = self.outputFluxes,
			objectiveValue = self.objectiveValue,
			objectiveComponents = self.objectiveComponents,
			)

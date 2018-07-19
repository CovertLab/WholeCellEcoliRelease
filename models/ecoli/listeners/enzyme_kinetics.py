#!/usr/bin/env python

"""
EnzymeKinetics

EnzymeKinetics listener. Tracks information about enzyme kinetics.

@organization: Covert Lab, Department of Bioengineering, Stanford University
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
	def initialize(self, sim, sim_data):
		super(EnzymeKinetics, self).initialize(sim, sim_data)

		self.metabolism = sim.processes["Metabolism"]
		self.metaboliteIDs = sorted(sim_data.process.metabolism.concDict)
		self.nConstrainedReactions = len(self.metabolism.kineticsConstrainedReactions)

		# Get metabolite names similar to how it's done in the metabolism process
		self.metaboliteNamesFromNutrients = set()
		for time, nutrientsLabel in sim_data.external_state.environment.nutrients_time_series[
			sim_data.external_state.environment.nutrients_time_series_label]:
				self.metaboliteNamesFromNutrients.update(
					sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
						nutrientsLabel, sim_data.process.metabolism.nutrientsToInternalConc
						)
					)
		self.metaboliteNamesFromNutrients = sorted(self.metaboliteNamesFromNutrients)

	# Allocate memory
	# In case things are of unknown size, write them here
	# Dummy values for what will be writen to output table
	# prep variables with zeros or NaNs, with correct size
	# to be filled later
	def allocate(self):
		super(EnzymeKinetics, self).allocate()
		self.metaboliteCountsInit = np.zeros(len(self.metaboliteNamesFromNutrients), np.float64)
		self.metaboliteCountsFinal = np.zeros(len(self.metaboliteNamesFromNutrients), np.float64)
		self.metaboliteConcentrations = np.zeros(len(self.metaboliteNamesFromNutrients), np.float64)
		self.enzymeIDs = self.metabolism.kineticsEnzymesList
		self.enzymeCountsInit = np.zeros(len(self.metabolism.kineticsEnzymesList), np.float64)
		self.countsToMolar = np.zeros(1, np.float64)
		self.targetFluxes = np.zeros(self.nConstrainedReactions, np.float64)
		self.actualFluxes = np.zeros(self.nConstrainedReactions, np.float64)
		self.reactionConstraint = np.zeros(self.nConstrainedReactions, np.int)

	def update(self):
		pass

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			enzymeIDs = self.enzymeIDs,
			metaboliteNames = self.metaboliteNamesFromNutrients,
			constrainedReactions = self.metabolism.kineticsConstrainedReactions,
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			metaboliteCountsInit = self.metaboliteCountsInit,
			metaboliteCountsFinal = self.metaboliteCountsFinal,
			metaboliteConcentrations = self.metaboliteConcentrations,
			countsToMolar = self.countsToMolar,
			enzymeCountsInit = self.enzymeCountsInit,
			targetFluxes = self.targetFluxes,
			actualFluxes = self.actualFluxes,
			reactionConstraint = self.reactionConstraint,
			)

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
	def initialize(self, sim, sim_data):
		super(FBAResults, self).initialize(sim, sim_data)

		self.metabolism = sim.processes["Metabolism"]

		self.objectiveValue = 0.0

		self.metaboliteNamesFromNutrients = set()
		for time, nutrientsLabel in sim_data.external_state.environment.nutrients_time_series[
			sim_data.external_state.environment.nutrients_time_series_label
			]:

			self.metaboliteNamesFromNutrients.update(
				sim_data.process.metabolism.concentrationUpdates.concentrationsBasedOnNutrients(
					nutrientsLabel, sim_data.process.metabolism.nutrientsToInternalConc
					)
				)
		self.metaboliteNamesFromNutrients = sorted(self.metaboliteNamesFromNutrients)


		if self.metabolism.run_flux_sensitivity:
			self.sensitivity_reactions = self.metabolism.kineticsConstrainedReactions + ['All enabled']
		else:
			self.sensitivity_reactions = []
		self.n_rxns_sensitivity = len(self.sensitivity_reactions)

	# Allocate memory
	def allocate(self):
		super(FBAResults, self).allocate()

		fba = self.metabolism.fba

		self.reactionIDs = fba.getReactionIDs()
		self.externalMoleculeIDs = fba.getExternalMoleculeIDs()
		self.outputMoleculeIDs = fba.getOutputMoleculeIDs()
		self.kineticTargetFluxNames = fba.getKineticTargetFluxNames()
		self.homeostaticTargetMolecules = fba.getHomeostaticTargetMolecules()

		self.reactionFluxes = np.zeros(len(self.reactionIDs), np.float64)
		self.externalExchangeFluxes = np.zeros(len(self.externalMoleculeIDs), np.float64)
		self.shadowPrices = np.zeros(len(self.outputMoleculeIDs), np.float64)
		self.reducedCosts = np.zeros(len(self.reactionIDs), np.float64)
		self.homeostaticObjectiveValues = np.zeros(len(self.homeostaticTargetMolecules))
		self.kineticObjectiveValues = np.zeros(len(self.kineticTargetFluxNames))
		self.deltaMetabolites = np.zeros(len(self.metaboliteNamesFromNutrients), np.float64)
		self.targetConcentrations = np.zeros(len(self.homeostaticTargetMolecules))

		self.succinate_flux_sensitivity = np.zeros(self.n_rxns_sensitivity)
		self.isocitrate_flux_sensitivity = np.zeros(self.n_rxns_sensitivity)


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			reactionIDs = list(self.reactionIDs),
			externalMoleculeIDs = self.externalMoleculeIDs,
			outputMoleculeIDs = self.outputMoleculeIDs,
			homeostaticTargetMolecules = self.homeostaticTargetMolecules,
			kineticTargetFluxNames = self.kineticTargetFluxNames,
			metaboliteNames = self.metaboliteNamesFromNutrients,
			sensitivity_reactions = self.sensitivity_reactions,
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			reactionFluxes = self.reactionFluxes,
			externalExchangeFluxes = self.externalExchangeFluxes,
			shadowPrices = self.shadowPrices,
			reducedCosts = self.reducedCosts,
			objectiveValue = self.objectiveValue,
			homeostaticObjectiveValues = self.homeostaticObjectiveValues,
			kineticObjectiveValues = self.kineticObjectiveValues,
			deltaMetabolites = self.deltaMetabolites,
			targetConcentrations = self.targetConcentrations,
			succinate_flux_sensitivity = self.succinate_flux_sensitivity,
			isocitrate_flux_sensitivity = self.isocitrate_flux_sensitivity,
			)

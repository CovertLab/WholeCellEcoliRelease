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

		# exchange with environment
		self.externalExchangeMolecules = self._external_states['Environment']._moleculeIDs

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			reactionIDs = list(self.reactionIDs),
			externalMoleculeIDs = self.externalMoleculeIDs,
			outputMoleculeIDs = self.outputMoleculeIDs,
			homeostaticTargetMolecules = self.homeostaticTargetMolecules,
			kineticTargetFluxNames = self.kineticTargetFluxNames,
			metaboliteNames = self.metaboliteNamesFromNutrients,
			externalExchangeMolecules = self.externalExchangeMolecules
			)


	def tableAppend(self, tableWriter):
		# save if importExchangeMolecules are constrained
		import_constraint = []
		for molecule_id in self.externalExchangeMolecules:
			if molecule_id in self.metabolism.exchange_data['importConstrainedExchangeMolecules'].keys():
				import_constraint.append(True)
			else:
				import_constraint.append(False)

		# save if externalExchangeMolecules are in importExchangeMolecules
		import_exchange = []
		for molecule_id in self.externalExchangeMolecules:
			if molecule_id in self.metabolism.exchange_data['importExchangeMolecules']:
				import_exchange.append(True)
			else:
				import_exchange.append(False)

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
			importConstraint = import_constraint,
			importExchange = import_exchange,
			)

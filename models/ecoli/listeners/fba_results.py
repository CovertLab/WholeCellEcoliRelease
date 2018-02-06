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

		self.countUnits = "counts"

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(FBAResults, self).initialize(sim, sim_data)

		self.metabolism = sim.processes["Metabolism"]

		self.objectiveValue = 0.0

		self.metaboliteNamesFromNutrients = set()
		for time, nutrientsLabel in sim_data.nutrientsTimeSeries[sim_data.nutrientsTimeSeriesLabel]:

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

		self.reactionIDs = fba.reactionIDs()
		self.externalMoleculeIDs = fba.externalMoleculeIDs()
		self.outputMoleculeIDs = fba.outputMoleculeIDs()
		self.kineticTargetFluxNames = fba.kineticTargetFluxNames()
		self.homeostaticTargetMolecules = fba.homeostaticTargetMolecules()

		self.deltaMetabolites = np.zeros(len(self.metaboliteNamesFromNutrients), np.float64)

		self.reactionFluxes = np.zeros(len(self.reactionIDs), np.float64)
		self.externalExchangeFluxes = np.zeros(len(self.externalMoleculeIDs), np.float64)
		self.outputFluxes = np.zeros(len(self.outputMoleculeIDs), np.float64)
		self.rowDualValues = np.zeros(len(self.outputMoleculeIDs), np.float64)
		self.columnDualValues = np.zeros(len(self.reactionIDs), np.float64)
		self.kineticObjectiveValues = np.zeros(len(self.kineticTargetFluxNames))
		self.homeostaticObjectiveValues = np.zeros(len(self.homeostaticTargetMolecules))
		self.homeostaticObjectiveWeight = np.ones(1, np.float64)
		self.targetConcentrations = np.zeros(len(self.homeostaticTargetMolecules))


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			reactionIDs = list(self.reactionIDs),
			externalMoleculeIDs = self.externalMoleculeIDs,
			outputMoleculeIDs = self.outputMoleculeIDs,
			homeostaticTargetMolecules = self.homeostaticTargetMolecules,
			kineticTargetFluxNames = self.kineticTargetFluxNames,
			metaboliteNames = self.metaboliteNamesFromNutrients,
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			reactionFluxes = self.reactionFluxes,
			externalExchangeFluxes = self.externalExchangeFluxes,
			outputFluxes = self.outputFluxes,
			rowDualValues = self.rowDualValues,
			columnDualValues = self.columnDualValues,
			objectiveValue = self.objectiveValue,
			kineticObjectiveValues = self.kineticObjectiveValues,
			homeostaticObjectiveValues = self.homeostaticObjectiveValues,
			homeostaticObjectiveWeight = self.homeostaticObjectiveWeight,
			deltaMetabolites = self.deltaMetabolites,
			targetConcentrations = self.targetConcentrations,
			)

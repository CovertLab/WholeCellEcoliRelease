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

		self.metabolism = sim.processes["Metabolism"].model

		self.objectiveValue = 0.0

		# exchange with environment
		self.all_external_exchange_molecules = sim_data.external_state.all_external_exchange_molecules

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
		self.deltaMetabolites = np.zeros(len(self.metabolism.metaboliteNamesFromNutrients), np.float64)
		self.targetConcentrations = np.zeros(len(self.homeostaticTargetMolecules))

		# exchange with environment
		self.import_constraint = [False] * len(self.all_external_exchange_molecules)
		self.import_exchange = [False] * len(self.all_external_exchange_molecules)


	def tableCreate(self, tableWriter):
		subcolumns = {
			'reactionFluxes': 'reactionIDs',
			'externalExchangeFluxes': 'externalMoleculeIDs',
			'shadowPrices': 'outputMoleculeIDs',
			'reducedCosts': 'reactionIDs',
			'homeostaticObjectiveValues': 'homeostaticTargetMolecules',
			'kineticObjectiveValues': 'kineticTargetFluxNames',
			'deltaMetabolites': 'metaboliteNames',
			'targetConcentrations': 'homeostaticTargetMolecules',
			'importConstraint': 'all_external_exchange_molecules',
			'importExchange': 'all_external_exchange_molecules'}

		tableWriter.writeAttributes(
			reactionIDs = list(self.reactionIDs),
			externalMoleculeIDs = self.externalMoleculeIDs,
			outputMoleculeIDs = self.outputMoleculeIDs,
			homeostaticTargetMolecules = self.homeostaticTargetMolecules,
			kineticTargetFluxNames = self.kineticTargetFluxNames,
			metaboliteNames = self.metabolism.metaboliteNamesFromNutrients,
			all_external_exchange_molecules = self.all_external_exchange_molecules,
			subcolumns = subcolumns)


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
			importConstraint = self.import_constraint,
			importExchange = self.import_exchange,
			)

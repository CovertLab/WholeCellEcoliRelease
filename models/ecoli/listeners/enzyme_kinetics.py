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
	def initialize(self, sim, sim_data):
		super(EnzymeKinetics, self).initialize(sim, sim_data)

		self.metabolism = sim.processes["Metabolism"]
#		self.reactionRateInfo = sim_data.process.metabolism.reactionRateInfo
		self.metaboliteIDs = sorted(sim_data.process.metabolism.concDict)
		self.nConstrainedReactions = len(sim_data.process.metabolism.constrainedReactionList)

		# Get metabolite names similar to how it's done in the metabolism process
		self.metaboliteNamesFromNutrients = set()
		for time, nutrientsLabel in sim_data.nutrientsTimeSeries[sim_data.nutrientsTimeSeriesLabel]:
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
#		self.baseRates = np.zeros(len(self.metabolism.fba.reactionIDs()), np.float64)
#		self.reactionKineticPredictions = np.zeros(len(self.metabolism.allRateReactions), np.float64)
#		self.allConstraintsLimits = np.zeros(len(self.reactionRateInfo), np.float64)
#		self.reactionIDs = self.metabolism.fba.reactionIDs()
#		self.kineticTargetFluxNames = self.metabolism.fba.kineticTargetFluxNames()
#		self.kineticOneSidedTargets = self.metabolism.fba.kineticOneSidedTargetFluxNames()
#		self.kineticTargetFluxes = np.zeros(len(self.metabolism.fba.kineticTargetFluxNames()), np.float64)
#		self.kineticTargetErrors = np.zeros(len(self.metabolism.fba.kineticTargetFluxNames()), np.float64)
#		self.kineticTargetRelativeDifferences = np.zeros(len(self.metabolism.fba.kineticTargetFluxNames()), np.float64)
#		self.overconstraintMultiples = np.zeros(len(self.metabolism.fba.reactionFluxes()[self.metabolism.allRateIndices]), np.float64)
#		self.constraintIDs = self.metabolism.constraintIDs
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
#			reactionIDs = list(self.reactionIDs),
#			constraintIDs = self.constraintIDs,
			enzymeIDs = self.enzymeIDs,
#			constraintToReactionDict = self.constraintToReactionDict,
#			kineticTargetFluxNames = self.kineticTargetFluxNames,
#			kineticOneSidedTargets = self.kineticOneSidedTargets,
			metaboliteNames = self.metaboliteNamesFromNutrients,
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
#			baseRates = self.baseRates,
#			reactionKineticPredictions = self.reactionKineticPredictions,
#			allConstraintsLimits = self.allConstraintsLimits,
#			kineticTargetFluxes = self.kineticTargetFluxes,
#			kineticTargetErrors = self.kineticTargetErrors,
#			kineticTargetRelativeDifferences = self.kineticTargetRelativeDifferences,
#			overconstraintMultiples = self.overconstraintMultiples,
			metaboliteCountsInit = self.metaboliteCountsInit,
			metaboliteCountsFinal = self.metaboliteCountsFinal,
			metaboliteConcentrations = self.metaboliteConcentrations,
			countsToMolar = self.countsToMolar,
			enzymeCountsInit = self.enzymeCountsInit,
			targetFluxes = self.targetFluxes,
			actualFluxes = self.actualFluxes,
			reactionConstraint = self.reactionConstraint,
			)

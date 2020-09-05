"""
FBAResults

Records dynamics of FBA output.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/24/2014
"""

from __future__ import absolute_import, division, print_function

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
		self.media_id = sim.external_states['Environment'].current_media_id

		self.objectiveValue = 0.0

		# exchange with environment
		self.all_external_exchange_molecules = sim_data.external_state.all_external_exchange_molecules
		self.conc_update_molecules = sim.processes["Metabolism"].conc_update_molecules

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

		# Args for metabolism functions
		self.conc_updates = np.zeros(len(self.conc_update_molecules))
		self.catalyst_counts = np.zeros(len(self.metabolism.catalyst_ids))
		self.translation_gtp = 0.
		self.coefficient = 0.
		self.unconstrained_molecules = [False] * len(self.all_external_exchange_molecules)
		self.constrained_molecules = [False] * len(self.all_external_exchange_molecules)
		self.uptake_constraints = [np.nan] * len(self.all_external_exchange_molecules)

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
			'importExchange': 'all_external_exchange_molecules',
			'conc_updates': 'conc_update_molecules',
			'catalyst_counts': 'catalyst_ids',
			}

		tableWriter.writeAttributes(
			reactionIDs=list(self.reactionIDs),
			externalMoleculeIDs=self.externalMoleculeIDs,
			outputMoleculeIDs=self.outputMoleculeIDs,
			homeostaticTargetMolecules=self.homeostaticTargetMolecules,
			kineticTargetFluxNames=self.kineticTargetFluxNames,
			metaboliteNames=self.metabolism.metaboliteNamesFromNutrients,
			all_external_exchange_molecules=self.all_external_exchange_molecules,
			conc_update_molecules=self.conc_update_molecules,
			catalyst_ids=self.metabolism.catalyst_ids,
			subcolumns=subcolumns,
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time=self.time(),
			simulationStep=self.simulationStep(),
			reactionFluxes=self.reactionFluxes,
			externalExchangeFluxes=self.externalExchangeFluxes,
			shadowPrices=self.shadowPrices,
			reducedCosts=self.reducedCosts,
			objectiveValue=self.objectiveValue,
			homeostaticObjectiveValues=self.homeostaticObjectiveValues,
			kineticObjectiveValues=self.kineticObjectiveValues,
			deltaMetabolites=self.deltaMetabolites,
			targetConcentrations=self.targetConcentrations,
			media_id=self.media_id,
			conc_updates=self.conc_updates,
			catalyst_counts=self.catalyst_counts,
			translation_gtp=self.translation_gtp,
			coefficient=self.coefficient,
			unconstrained_molecules=self.unconstrained_molecules,
			constrained_molecules=self.constrained_molecules,
			uptake_constraints=self.uptake_constraints,
			)

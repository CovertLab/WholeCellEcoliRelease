"""
EnzymeKinetics

EnzymeKinetics listener. Tracks information about enzyme kinetics.

"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener


class EnzymeKinetics(wholecell.listeners.listener.Listener):
	""" EnzymeKinetics """

	_name = 'EnzymeKinetics'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(EnzymeKinetics, self).__init__(*args, **kwargs)

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(EnzymeKinetics, self).initialize(sim, sim_data)

		self.metabolism = sim.processes["Metabolism"].model
		self.n_constrained_reactions = len(self.metabolism.kinetics_constrained_reactions)
		self.n_metabolites = len(self.metabolism.metaboliteNamesFromNutrients)

		self.constraint_is_kcat_only = sim_data.process.metabolism.constraint_is_kcat_only[
			self.metabolism.active_constraints_mask].tolist()

	# Allocate memory
	# In case things are of unknown size, write them here
	# Dummy values for what will be writen to output table
	# prep variables with zeros or NaNs, with correct size
	# to be filled later
	def allocate(self):
		super(EnzymeKinetics, self).allocate()
		self.metaboliteCountsInit = np.zeros(self.n_metabolites, np.float64)
		self.metaboliteCountsFinal = np.zeros(self.n_metabolites, np.float64)
		self.enzymeIDs = self.metabolism.kinetic_constraint_enzymes
		self.enzymeCountsInit = np.zeros(len(self.metabolism.kinetic_constraint_enzymes), np.float64)
		self.countsToMolar = np.zeros(1, np.float64)
		self.targetFluxes = np.zeros(self.n_constrained_reactions, np.float64)
		self.targetFluxesUpper = np.zeros(self.n_constrained_reactions, np.float64)
		self.targetFluxesLower = np.zeros(self.n_constrained_reactions, np.float64)
		self.actualFluxes = np.zeros(self.n_constrained_reactions, np.float64)

	def update(self):
		pass

	def tableCreate(self, tableWriter):
		subcolumns = {
			'metaboliteCountsInit': 'metaboliteNames',
			'metaboliteCountsFinal': 'metaboliteNames',
			'enzymeCountsInit': 'enzymeIDs',
			'targetFluxes': 'constrainedReactions',
                        'targetFluxesUpper': 'constrainedReactions',
                        'targetFluxesLower': 'constrainedReactions',
			'actualFluxes': 'constrainedReactions'}

		tableWriter.writeAttributes(
			enzymeIDs = self.enzymeIDs,
			metaboliteNames = self.metabolism.metaboliteNamesFromNutrients,
			constrainedReactions = self.metabolism.kinetics_constrained_reactions,
			kineticsConstrainedReactions = self.metabolism.kinetics_constrained_reactions,
			constraint_is_kcat_only = self.constraint_is_kcat_only,
			subcolumns = subcolumns)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			metaboliteCountsInit = self.metaboliteCountsInit,
			metaboliteCountsFinal = self.metaboliteCountsFinal,
			countsToMolar = self.countsToMolar,
			enzymeCountsInit = self.enzymeCountsInit,
			targetFluxes = self.targetFluxes,
                        targetFluxesUpper = self.targetFluxesUpper,
                        targetFluxesLower = self.targetFluxesLower,
			actualFluxes = self.actualFluxes,
			)

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
		self.n_constrained_reactions = len(self.metabolism.kinetics_constrained_reactions)

		# flux targets from boundary
		self.n_boundary_constrained_reactions = len(self.metabolism.boundary_constrained_reactions)
		self.n_all_constrained_reactions = self.n_constrained_reactions + self.n_boundary_constrained_reactions
		self.n_metabolites = len(self.metabolism.metaboliteNamesFromNutrients)

	# Allocate memory
	# In case things are of unknown size, write them here
	# Dummy values for what will be writen to output table
	# prep variables with zeros or NaNs, with correct size
	# to be filled later
	def allocate(self):
		super(EnzymeKinetics, self).allocate()
		self.metaboliteCountsInit = np.zeros(self.n_metabolites, np.float64)
		self.metaboliteCountsFinal = np.zeros(self.n_metabolites, np.float64)
		self.metaboliteConcentrations = np.zeros(self.n_metabolites, np.float64)
		self.enzymeIDs = self.metabolism.kineticsEnzymesList
		self.enzymeCountsInit = np.zeros(len(self.metabolism.kineticsEnzymesList), np.float64)
		self.countsToMolar = np.zeros(1, np.float64)
		self.targetFluxes = np.zeros(self.n_all_constrained_reactions, np.float64)
		self.actualFluxes = np.zeros(self.n_all_constrained_reactions, np.float64)

		# reactionConstraint is only for kinetic constrained reactions, without boundary constrained reactions
		self.reactionConstraint = np.zeros(self.n_constrained_reactions, np.int)

	def update(self):
		pass

	def tableCreate(self, tableWriter):
		subcolumns = {
			'metaboliteCountsInit': 'metaboliteNames',
			'metaboliteCountsFinal': 'metaboliteNames',
			'metaboliteConcentrations': 'metaboliteNames',
			'enzymeCountsInit': 'enzymeIDs',
			'targetFluxes': 'constrainedReactions',
			'actualFluxes': 'constrainedReactions'}

		tableWriter.writeAttributes(
			enzymeIDs = self.enzymeIDs,
			metaboliteNames = self.metabolism.metaboliteNamesFromNutrients,
			constrainedReactions = self.metabolism.all_constrained_reactions,
			kineticsConstrainedReactions = self.metabolism.kinetics_constrained_reactions,
			boundaryConstrainedReactions = self.metabolism.boundary_constrained_reactions,
			subcolumns = subcolumns)


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

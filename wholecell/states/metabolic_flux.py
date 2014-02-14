#!/usr/bin/env python

"""
MetabolicFlux

MetabolicFlux state variable. Represents the instantaneous growth rate (fg/h) and metabolic reaction fluxes (reactions/s)

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

import numpy

import wholecell.states.state

class MetabolicFlux(wholecell.states.state.State):
	""" MetabolicFlux """

	# Constructor
	def __init__(self, *args, **kwargs):
		self.meta = {
			"id": "MetabolicFlux",
			"name": "MetabolicFlux",
			"dynamics": ["growth", "fluxes"],
			"units": {
				"growth": "fg/h",
				"fluxes": "reactions/s"
				}
		}
		# State, process references
		self.moleculeCounts = None
		self.metabolism = None

		self.growth = None
		self.fluxes = None

		super(MetabolicFlux, self).__init__(*args, **kwargs)

	# Construct state-process graph
	def initialize(self, sim, kb):
		super(MetabolicFlux, self).initialize(sim, kb)

		self.moleculeCounts = sim.getState("MoleculeCounts")
		self.metabolism = sim.getProcess("Metabolism")

		self.reactionIds = [reaction["id"] for reaction in kb.reactions]
		self.reactionNames = [reaction["name"] for reaction in kb.reactions]

	# Allocate memory
	def allocate(self):
		super(MetabolicFlux, self).allocate()

		self.growth = numpy.zeros(1)
		self.fluxes = numpy.zeros(len(self.reactionIds)) # Note: Probably not the right size

	# Calculate initial conditions
	def calcInitialConditions(self):
		mc = self.moleculeCounts
		met = self.metabolism

		# TODO: Fix
		# bounds = met.calcFluxBounds(
		# 	mc.counts[numpy.unravel_index(met.metabolite.mapping, mc.counts.shape)], mc.counts[numpy.unravel_index(met.enzyme.mapping, mc.counts.shape)]
		# 	)
		# bounds = met.calcFluxBounds(
		# 	mc.counts[numpy.unravel_index(met.metabolite.mapping, mc.counts.shape)], numpy.zeros(1), applyKineticBounds = False, applyExchangeBounds = False
		# 	)

		# self.growth, self.fluxes = met.calcGrowthRate(bounds)[0:2]

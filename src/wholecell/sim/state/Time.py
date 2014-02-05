#!/usr/bin/env python

"""
Time

Time state variable. Represents the current time lapsed since the start of the simulation in seconds.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

import numpy

import wholecell.sim.state.State

class Time(wholecell.sim.state.State.State):
	""" Time """

	# Constructor
	def __init__(self, *args, **kwargs):
		self.meta = {
			"id": "Time",
			"name": "Time",
			"dynamics": ["value"],
			"units": {"value": "s"}
		}

		self.value = None			# Simulation time in seconds

		super(Time, self).__init__(*args, **kwargs)


	def initialize(self, sim, kb):
		super(Time, self).initialize(sim, kb)

		self.sim = sim # Time needs a reference to the simulation to determine what step the simulation is on
		self.timeStepSec = sim.timeStepSec


	# Allocate memory
	def allocate(self):
		super(Time, self).allocate()

		self.value = numpy.zeros(1)


	def calculate(self):
		self.value = self.timeStepSec * (self.sim.initialStep + self.sim.simulationStep)

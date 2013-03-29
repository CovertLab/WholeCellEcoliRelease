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

	meta = {
		"id": "Time",
		"name": "Time",
		"dynamics": ["value"],
		"units": {"value": "s"}
	}

	def __init__(self, *args, **kwargs):
		self.value = None
		super(Time, self).__init__(*args, **kwargs)

	def allocate(self):
		super(Time, self).allocate()

		self.value = numpy.zeros(1)

	def calcInitialConditions(self):
		self.value = 0
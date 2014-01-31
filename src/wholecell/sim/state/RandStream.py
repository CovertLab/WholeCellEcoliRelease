#!/usr/bin/env python

"""
RandStream

Pseudorandom number generator state.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/29/2013
"""

import numpy

import wholecell.sim.state.State

class RandStream(wholecell.sim.state.State.State):
	""" RandStream """

	# Constructor
	def __init__(self, *args, **kwargs):
		self.meta = {
			"id": "RandStream",
			"name": "RandStream",
			"dynamics": ["value"],
			# "units": {"value": "s"}
		}

		self.value = None
		super(RandStream, self).__init__(*args, **kwargs)

	
	def calculate(self):
		self.value = self.randStream.state

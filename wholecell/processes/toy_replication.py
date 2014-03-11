#!/usr/bin/env python

"""
ToyReplication

A toy process for developing interactions between processes and the Chromosome
State.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 3/11/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class ToyReplication(wholecell.processes.process.Process):
	""" ToyReplication """

	# Constructor
	def __init__(self):
		self.meta = {
			"id": "ToyReplication",
			"name": "ToyReplication",
		}

		super(ToyReplication, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(ToyReplication, self).initialize(sim, kb)

		pass


	def calculateRequest(self):
		pass


	# Calculate temporal evolution
	def evolveState(self):
		pass

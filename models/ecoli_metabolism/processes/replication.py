#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process

class Replication(wholecell.processes.process.Process):
	""" Replication """

	_name = "Replication"

	# Constructor
	def __init__(self):
		super(Replication, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Replication, self).initialize(sim, kb)

		# Load constants

		# Create views on state


	def calculateRequest(self):
		pass


	def evolveState(self):
		pass

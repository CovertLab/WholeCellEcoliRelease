#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process

class Maintenance(wholecell.processes.process.Process):
	""" Maintenance """

	_name = "Maintenance"

	# Constructor
	def __init__(self):
		super(Maintenance, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Maintenance, self).initialize(sim, kb)

		# Load constants

		# Create views on state


	def calculateRequest(self):
		pass


	def evolveState(self):
		pass

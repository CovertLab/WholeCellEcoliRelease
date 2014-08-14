#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process

class Metabolism(wholecell.processes.process.Process):
	""" Metabolism """

	_name = "Metabolism"

	# Constructor
	def __init__(self):
		super(Metabolism, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Metabolism, self).initialize(sim, kb)

		# Load constants

		# Create views on state


	def calculateRequest(self):
		pass


	def evolveState(self):
		pass

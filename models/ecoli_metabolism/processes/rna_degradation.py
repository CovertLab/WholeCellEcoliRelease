#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process

class RnaDegradation(wholecell.processes.process.Process):
	""" RnaDegradation """

	_name = "RnaDegradation"

	# Constructor
	def __init__(self):
		super(RnaDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(RnaDegradation, self).initialize(sim, kb)

		# Load constants

		# Create views on state


	def calculateRequest(self):
		pass


	def evolveState(self):
		pass

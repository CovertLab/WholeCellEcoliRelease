#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process

class ProteinDegradation(wholecell.processes.process.Process):
	""" ProteinDegradation """

	_name = "ProteinDegradation"

	# Constructor
	def __init__(self):
		super(ProteinDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(ProteinDegradation, self).initialize(sim, kb)

		# Load constants

		# Create views on state


	def calculateRequest(self):
		pass


	def evolveState(self):
		pass
